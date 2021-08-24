/*
 *  inhomogeneous_poisson_generator.cpp
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

#include "stream_inhomogeneous_poisson_generator.h"

// C++ includes:

// Includes from libnestutil:

// Includes from nestkernel:
#include "event_delivery_manager_impl.h"
#include "kernel_manager.h"

// Includes from sli:

/* ----------------------------------------------------------------
 * Update function and event hook
 * ---------------------------------------------------------------- */
void
nest::stream_inhomogeneous_poisson_generator::update( Time const& origin, const long from, const long to )
{
  // not update by nest itself only by the background
}


void
nest::stream_inhomogeneous_poisson_generator::event_hook( DSSpikeEvent& e )
{
  poisson_distribution::param_type param( rate_ * V_.h_ );
  long n_spikes = V_.poisson_dist_( get_vp_specific_rng( get_thread() ), param );

  if ( n_spikes > 0 ) // we must not send events with multiplicity 0
  {
    e.set_multiplicity( n_spikes );
    e.get_receiver().handle( e );
  }
}

/* ----------------------------------------------------------------
 * Other functions
 * ---------------------------------------------------------------- */
void
nest::stream_inhomogeneous_poisson_generator::set_data_from_stream_stimulating_backend( double& rate )
{
  rate_ = rate;
  // create spikes
  if ( rate_ > 0 )
  {
    DSSpikeEvent se;
    kernel().event_delivery_manager.send( *this, se, 0.0 );
  }
}