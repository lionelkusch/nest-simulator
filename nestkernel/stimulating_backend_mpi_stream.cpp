/*
 *  stimulating_backend_mpi.cpp
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
#include <fstream>

// Includes from nestkernel:
#include "stimulating_backend_mpi_stream.h"
#include "kernel_manager.h"
#include "stimulating_device.h"


void
nest::StimulatingBackendMPIStream::prepare()
{
  if ( not enrolled_ )
  {
    return;
  }

  if ( prepared_ )
  {
    throw BackendPrepared( "StimulatingBackendMPI" );
  }

  // need to be run only by the master thread : it is the case because this part is not running in parallel
  thread thread_id_master = kernel().vp_manager.get_thread_id();
  // Create the connection with MPI
  // 1) take all the ports of the connections
  // get port and update the list of device only for master
  for ( auto& it_device : devices_[ thread_id_master ] )
  {
    // add the link between MPI communicator and the device (devices can share the same MPI communicator)
    std::string port_name;
    get_port( it_device.second.second, &port_name );
    auto comm_it = commMap_.find( port_name );
    MPI_Comm* comm;
    if ( comm_it != commMap_.end() )
    {
      // it's not a new communicator
      comm = std::get< 0 >( comm_it->second );
      // add the id of the device if there are a connection with the device.
      std::get< 1 >( comm_it->second )->push_back( it_device.second.second->get_node_id() );
      std::get< 2 >( comm_it->second )[ thread_id_master ] += 1;
    }
    else
    {
      // create a new MPI communicator to communicate with the external MPI process.
      // Only the master thread uses the MPI functions of this new communicator.
      // This is because the management of threads here is using MPI_THREAD_FUNNELED (see mpi_manager.cpp:119).
      comm = new MPI_Comm;
      auto vector_id_device = new std::vector< int >; // vector of ID device for the rank
      int* vector_nb_device_th{ new int[ kernel().vp_manager.get_num_threads() ]{} }; // number of device by thread
      std::fill_n( vector_nb_device_th, kernel().vp_manager.get_num_threads(), 0 );
      // add the id of the device if there is a connection with the device.
      if ( kernel().connection_manager.get_device_connected(
        thread_id_master, it_device.second.second->get_local_device_id() ) )
      {
        vector_id_device->push_back( it_device.second.second->get_node_id() );
        vector_nb_device_th[ thread_id_master ] += 1;
      }
      std::tuple< MPI_Comm*, std::vector< int >*, int* > comm_count =
        std::make_tuple( comm, vector_id_device, vector_nb_device_th );
      commMap_.insert( std::make_pair( port_name, comm_count ) );
    }
    it_device.second.first = comm;
  }

  // 2) connect the master thread to the MPI process it needs to be connected to
  for ( auto& it_comm : commMap_ )
  {
    MPI_Comm_connect( it_comm.first.data(),
                      MPI_INFO_NULL,
                      0,
                      MPI_COMM_WORLD,
                      std::get< 0 >( it_comm.second ) ); // should use the status for handle error
    std::ostringstream msg;
    msg << "Connect to " << it_comm.first.data() << "\n";
    LOG( M_INFO, "MPI Input connect", msg.str() );
  }
  step_ = 1;
  double** data = { new double*[ commMap_.size() ]{} };
  data_ = &data;
}

void
nest::StimulatingBackendMPIStream::pre_run_hook()
{
  // nothing to do
}

void
nest::StimulatingBackendMPIStream::pre_step_hook()
{
#pragma omp master
  {
    int index = 0;
    // receive all the information from all the MPI connections
    for ( auto& it_comm : commMap_ )
    {
      // Receive the new value of device
      MPI_Status status_mpi;
      int shape = { *std::get<2>(it_comm.second) };
      double* data_receive{ new double[ shape ]{} };
      MPI_Recv( data_receive, shape, MPI_DOUBLE, 0, step_, *std::get< 0 >( it_comm.second), &status_mpi );
      (*data_)[index] = data_receive;
      index += 1;
    }
    step_+=1;
  }
#pragma omp barrier
  // Each thread updates its own devices.
  int index_it = 0;
  for ( auto& it_comm : commMap_ )
  {
    update_device( *std::get< 1 >( it_comm.second ), ( *data_)[ index_it ] );
    index_it += 1;
  }
#pragma omp barrier
#pragma omp master
  {
   for (int i=0; i < (int) commMap_.size(); i++){
     delete[] (*data_)[i];
   }
  }
}

void
nest::StimulatingBackendMPIStream::post_run_hook()
{
  // nothing to do
}

void
nest::StimulatingBackendMPIStream::cleanup()
{
// Disconnect all the MPI connection and send information about this disconnection
// Clean all the elements in the map
// disconnect MPI message
#pragma omp master
  {
    for ( auto& it_comm : commMap_ )
    {
      bool value[ 1 ] = { true };
      MPI_Send( value, 1, MPI_CXX_BOOL, 0, 2, *std::get< 0 >( it_comm.second ) );
      MPI_Barrier(*std::get< 0 >( it_comm.second ));
      MPI_Comm_disconnect( std::get< 0 >( it_comm.second ) );
      delete std::get< 0 >( it_comm.second );
      delete std::get< 1 >( it_comm.second );
      delete[] std::get< 2 >( it_comm.second );
      std::get< 2 >( it_comm.second ) = nullptr;
    }
    // clear map of devices
    commMap_.clear();
    thread thread_id_master = kernel().vp_manager.get_thread_id();
    for ( auto& it_device : devices_[ thread_id_master ] )
    {
      it_device.second.first = nullptr;
    }
    delete[] *data_;
  }
#pragma omp barrier
}


void
nest::StimulatingBackendMPIStream::update_device( std::vector< int >& devices_id, double* data )
{
  // if there are some data
  thread thread_id = kernel().vp_manager.get_thread_id();
  for ( int i = 0; i < ( int ) devices_id.size(); i++ )
  {
    devices_[ thread_id ].find( devices_id[ i ] )->second.second->set_data_from_stream_stimulating_backend( data[ i ] );
  }
}