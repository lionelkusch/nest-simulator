/*
 *  stimulating_backend_mpi.h
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

#ifndef STIMULATING_BACKEND_MPI_STREAM_H
#define STIMULATING_BACKEND_MPI_STREAM_H

#include "stimulating_backend_mpi.h"
#include "nest_types.h"
#include "nest_time.h"
#include <set>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <mpi.h>

/* BeginUserDocs: stimulating backend

.. _stimulating_backend_mpi:

Collect data from MPI communication for updating device
#######################################################

.. admonition:: Availability

   This stimulating backend is only available if NEST was compiled with
   :ref:`support for MPI <compile-with-mpi>`.

The `mpi` stimulating backend collects data from MPI channels and
updates stimulating device just before each run. This is useful for
co-simulation or for receiving data from an external software.

The creation of this MPI communication is based on the functions
'MPI_Comm_connect' and 'MPI_Comm_disconnect'. The port name is
read from a file for each device with this backend.
The name of the file needs to be specified according to the following
pattern:

::

   {data_path}/{data_prefix}{label}/{id_device}.txt

The ``data_path`` and ``data_prefix`` are global kernel properties,
while `label` is a property of the device in question and `id_device`
its node ID.
This path can only be set outside of a `Run` contexts (i.e.
after ``Prepare()` has been called, but ``Cleanup()`` has not).


Communication Protocol
++++++++++++++++++++++
The following protocol is used to exchange information between
both MPI processes. The protocol is described using the
following format for the MPI messages:
(value, number, type, source/destination, tag)

1) ``Prepare``  : Connection of MPI port include in one file (see below)
2) ``Run`` begin: Send start run (true, 1, CXX_BOOL, 0, 0)
3) ``Run`` begin: Send the id of the device to update (id_device, 1, INT, 0, 0)
3) ``Run`` begin: Receive shape of the data (shape, 1, INT, 0, 0)
4) ``Run`` begin: Receive the data for updating the device (data, shape, DOUBLE, 0, 0)
5) ``Run`` end  : Send at each ending of the run (true, 1, CXX_BOOL, 0, 1)
6) ``Cleanup``  : Send at this en of the simulation (true, 1, CXX_BOOL, 0, 2)

Data format
+++++++++++

The format of the data is depend of the stimulating device.

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
 * A simple input backend MPI implementation
 */
class StimulatingBackendMPIStream : public StimulatingBackendMPI
{
public:
  void prepare() override;

  void pre_run_hook() override;

  void post_run_hook() override;

  void pre_step_hook() override;

  void cleanup() override;

private:
  /**
   * Update all the devices with the data received
   * @param array_index : number of devices by thread
   * @param devices_id : the devices' ID ordered by thread
   * @param data : the data received for updating all the devices
   */
  void update_device( std::vector< int >& devices_id, double* data );

  int step_;
  double*** data_;
};

} // namespace

#endif // STIMULATING_BACKEND_MPI_STREAM_H
