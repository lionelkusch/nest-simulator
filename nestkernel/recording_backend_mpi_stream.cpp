/*
 *  recording_backend_mpi.cpp
 *
 *  This file is part of NEST.
 *
 *  Copyright (C) 2004 The NEST Initiative
 *
 *  NEST is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  NEST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with NEST.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

// C++ includes:
#include <iostream>


// Includes from nestkernel:
#include "recording_device.h"
#include "recording_backend_mpi.h"
#include "recording_backend_mpi_stream.h"
#include "exceptions.h"

void
nest::RecordingBackendMPIStream::initialize()
{
  auto nthreads = kernel().vp_manager.get_num_threads();
  std::vector< std::vector< int > > empty_vector( nthreads );
  buffer_stream_.swap( empty_vector );
  device_map devices( nthreads );
  devices_.swap( devices );
}

void
nest::RecordingBackendMPIStream::prepare()
{
  if ( not enrolled_ )
  {
    return;
  }

  if ( prepared_ )
  {
    throw BackendPrepared( "RecordingBackendMPI" );
  }
  prepared_ = true;
  thread thread_id_master = 0;
#pragma omp parallel default( none ) shared( thread_id_master )
  {
#pragma omp master
    {
      // Create the connection with MPI
      // 1) take all the ports of the connections
      // get port and update the list of devices
      thread_id_master = kernel().vp_manager.get_thread_id();
    }
  }
  int count_max = 0;
  for ( auto& it_device : devices_[ thread_id_master ] )
  {
    // add the link between MPI communicator and the device (devices can share the same MPI communicator)
    std::string port_name;
    get_port( std::get< 2 >( it_device.second ), &port_name );
    auto comm_it = commMap_.find( port_name );
    MPI_Comm* comm;
    int index_mpi;
    if ( comm_it != commMap_.end() )
    {
      comm = std::get< 1 >( comm_it->second );
      std::get< 2 >( comm_it->second ) += 1;
      index_mpi = std::get< 0 >( comm_it->second );
    }
    else
    {
      comm = new MPI_Comm;
      std::tuple< int, MPI_Comm*, int > comm_count = std::make_tuple( count_max, comm, 1 );
      commMap_.insert( std::make_pair( port_name, comm_count ) );
      index_mpi = count_max;
      count_max += 1;
    }
    std::get< 0 >( it_device.second ) = index_mpi;
    std::get< 1 >( it_device.second ) = comm;
  }

  // initialize the buffer
  for ( auto& thread_data : buffer_stream_ )
  {
    std::vector< int > data_comm( count_max );
    thread_data.swap( data_comm );
  }

  // 2) connect the thread to the MPI process it needs to be connected to
  for ( auto& it_comm : commMap_ )
  {
    MPI_Comm_connect( it_comm.first.data(),
      MPI_INFO_NULL,
      0,
      MPI_COMM_WORLD,
      std::get< 1 >( it_comm.second ) ); // should use the status for handle error
    std::ostringstream msg;
    msg << "Connect to " << it_comm.first.data() << "\n";
    LOG( M_INFO, "MPI Record connect", msg.str() );
  }
#pragma omp parallel default( none ) shared( thread_id_master )
  {
    // Update all the threads
    thread thread_id = kernel().vp_manager.get_thread_id();
    if ( thread_id != thread_id_master )
    {
      for ( auto& it_device : devices_[ thread_id ] )
      {
        auto device_data = devices_[ thread_id_master ].find( it_device.first );
        std::get< 0 >( it_device.second ) = std::get< 0 >( device_data->second );
        std::get< 1 >( it_device.second ) = std::get< 1 >( device_data->second );
      }
    }
  }
  step_=0;
}


void
nest::RecordingBackendMPIStream::pre_run_hook()
{
  // nothing to do
}


void
nest::RecordingBackendMPIStream::post_step_hook()
{
#pragma omp master
  {
    step_+=1;
    // Receive information of MPI process
    for ( auto& it_comm : commMap_ )
    {
      std::vector< int > data;
      for ( auto& data_thread : buffer_stream_ )
      {
        data.push_back(0);
        for ( auto& data_sample : data_thread )
        {
          data.back() += data_sample;
        }
      }
      send_data( std::get< 1 >( it_comm.second ),  data.data(), data.size(), step_);
    }
    // clear the buffer
    for ( auto& data_thread : buffer_stream_ )
    {
      for ( auto& data_sample : data_thread )
      {
        data_sample=0;
      }
    }
    // Send information about the end of the running part
  }
#pragma omp barrier
}

void
nest::RecordingBackendMPIStream::post_run_hook()
{
  // nothing to do
}

void
nest::RecordingBackendMPIStream::cleanup()
{
// Disconnect all the MPI connections and send information about this disconnection
// Clean all the elements in the map
// disconnect MPI
#pragma omp master
  {
    for ( auto& it_comm : commMap_ )
    {
      int  data[ std::get< 0 >(it_comm.second) ];
      send_data(std::get< 1 >( it_comm.second ),data,std::get< 0 >(it_comm.second),0);

      MPI_Barrier(*std::get< 1 >( it_comm.second ));
      MPI_Comm_disconnect( std::get< 1 >( it_comm.second ) );
      delete ( std::get< 1 >( it_comm.second ) );
    }
    // clear the buffer
    for ( auto& data_thread : buffer_stream_ )
    {
      data_thread.clear();
    }
    // clear map of device
    commMap_.clear();
    thread thread_id_master = kernel().vp_manager.get_thread_id();
    for ( auto& it_device : devices_[ thread_id_master ] )
    {
      std::get< 0 >( it_device.second ) = -1;
      std::get< 1 >( it_device.second ) = nullptr;
    }
  }
#pragma omp barrier
}

void
nest::RecordingBackendMPIStream::write( const RecordingDevice& device,
  const Event& event,
  const std::vector< double >&,
  const std::vector< long >& )
{
  // For each event send a message through the right MPI communicator
  const thread thread_id = kernel().get_kernel_manager().vp_manager.get_thread_id();
  const index recorder = device.get_node_id();

  auto it_devices = devices_[ thread_id ].find( recorder );
  if ( it_devices != devices_[ thread_id ].end() )
  {
    buffer_stream_[ thread_id ][ std::get< 0 >( it_devices->second ) ] += \
      reinterpret_cast < const SpikeEvent*> (&event)->get_multiplicity();
  }
  else
  {
    throw BackendPrepared( " Internal error " );
  }
}

/* ----------------------------------------------------------------
 * Parameter extraction and manipulation functions
 * ---------------------------------------------------------------- */
void
nest::RecordingBackendMPIStream::send_data( const MPI_Comm* comm, const int data[], const int size, int tag )
{
  // Send the size of data
  int shape = { size };
  // Receive the data ( for the moment only spike time )
  MPI_Send( data, shape, MPI_INT, 0, tag, *comm );
//  std::cout<<"NEST :"<<data[0]<<std::endl; std::cout.flush();
}
