/*
 *  recording_backend_mpi_stream.h
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

#ifndef RECORDING_BACKEND_MPI_STREAM_H
#define RECORDING_BACKEND_MPI_STREAM_H

#include "recording_backend_mpi.h"

/* BeginUserDocs: recording backend

.. _recording_backend_mpi_stream:

Send data with MPI
##################

.. admonition:: Availability

   This stimulating backend is only available if NEST was compiled with
   :ref:`support for MPI <compile-with-mpi>`.

The `mpi` recording backend sends collected data to a remote process
using MPI.

This backend will create a new MPI communicator (different from
MPI_Comm_World, which is used by NEST itself).  The creation of the
MPI communication is based on the functions 'MPI_Comm_connect' and
'MPI_Comm_disconnect'. The port name is read from a file for each
device with this backend. The file needs to be named according to the
following pattern:

::

   {data_path}/{data_prefix}{label}/{id_device}.txt

The ``data_path`` and ``data_prefix`` are global kernel properties,
while `label` is a property of the device in question and `id_device`
its node ID.
This path can only be set outside of a `Run` contexts (i.e.
after ``Prepare()` has been called, but ``Cleanup()`` has not).

Communication Protocol:
+++++++++++++++++++++++
The following protocol is used to exchange information between
both MPI processes. The protocol is described using the
following format for the MPI messages:
(value, number, type, source/destination, tag)

1) ``Prepare``  : Connection of MPI port included in the port_file (see below)
2) ``Run`` begin: Send at each beginning of the run (true, 1, CXX_BOOL, 0, 0)
3) ``Run`` end  : Receive at each ending of the run (true, 1, CXX_BOOL, 0, 0)
4) ``Run`` end  : Send shape of the data of the run (shape, 1,INT, 0, 0)
5) ``Run`` end  : Send data of the data of the run (data, shape, DOUBLE, 0, 0)
6) ``Run`` end  : Send at each ending of the run (true, 1, CXX_BOOL, 0, 1)
7) ``Cleanup``  : Send at this en of the simulation (true, 1, CXX_BOOL, 0, 2)

Data format
+++++++++++

The format of the sending data is an array of (id device, id node, time is ms).

Parameter summary
+++++++++++++++++

.. glossary::

 label
   Shared file with ports with the format path+label+id+.txt
   Where path = kernel().io_manager.get_data_path() and the
   label and id are the ones provided for the device.
EndUserDocs */

namespace nest
{

/**
 * A recording backend for sending information with MPI.
 */
class RecordingBackendMPIStream : public RecordingBackendMPI
{
public:
  void initialize() override;

  void prepare() override;

  void write( const RecordingDevice&, const Event&, const std::vector< double >&, const std::vector< long >& ) override;

  void pre_run_hook() override;

  void post_run_hook() override;

  void post_step_hook() override;

  void cleanup() override;

private:

  static void send_data( const MPI_Comm* comm, const int data[], int size, int tag );
    /**
   * Buffer for saving events before they are sent.
   * The buffer has 3 dimensions : thread_id, MPI_communicator_index and number of events
   * elements. The events elements are described as an array with three components: id of device, id of neurons and data
   * ( one double )
   */
  std::vector< std::vector < std::vector< int > > > buffer_stream_;

  int step_;
};

} // namespace

#endif // RECORDING_BACKEND_MPI_STREAM_H
