/*
 *  inhomogeneous_poisson_generator.h
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

#ifndef INHOMOGENEOUS_POISSON_GENERATOR_STREAM_H
#define INHOMOGENEOUS_POISSON_GENERATOR_STREAM_H

// C++ includes:
#include <vector>

// Includes from nestkernel:
#include "connection.h"
#include "device_node.h"
#include "event.h"
#include "nest.h"
#include "random_generators.h"
#include "ring_buffer.h"
#include "inhomogeneous_poisson_generator.h"

namespace nest
{
class stream_inhomogeneous_poisson_generator : public inhomogeneous_poisson_generator
{

public:
  void set_data_from_stream_stimulating_backend( double& rate ) override;
  StimulatingDevice::Type get_type() const override;

private:
  void event_hook( DSSpikeEvent& ) override;
  double rate_ = 0; //!< current amplitude
};

inline StimulatingDevice::Type
stream_inhomogeneous_poisson_generator::get_type() const
{
  return StimulatingDevice::Type::STREAM_SPIKE_GENERATOR;
}

} // namespace

#endif // INHOMOGENEOUS_POISSON_GENERATOR_STREAM_H
