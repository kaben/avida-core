/*
 *  unittests/core/Strand.cc
 *  avida-core
 *
 *  Created by Kaben on 7/31/13.
 *  Copyright 2013 Michigan State University. All rights reserved.
 *  http://avida.devosoft.org/
 *
 *
 *  This file is part of Avida.
 *
 *  Avida is free software; you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
 *  as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
 *
 *  Avida is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public License for more details.
 *
 *  You should have received a copy of the GNU Lesser General Public License along with Avida.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 *  Authors: Kaben G. Nanlohy <kaben.nanlohy@gmail.com>
 *
 */

#include "cpu/cStrand.h"

#include "gtest/gtest.h"


namespace nStrandBrainstorming {
  TEST(cStrand, strand_default_instantiation) {
    /*
     A default strand is essentially an empty string. Elsewhere in Avida,
     the convention is to initialize these default sequences to have length 1,
     so we're following that convention.
     */
    cStrand *strand = new cStrand();
    EXPECT_TRUE(strand);
    EXPECT_EQ(strand->GetSize(), 1);
    delete strand;
  }
  TEST(cStrand, strand_string_instantiation) {
    /*
     We're only going to test one other kind of instantiation:
     with a string argument. Since cStrand calls its superclass' constructors,
     which should be tested elsewhere,
     we're going to skip the rest of those tests here.
     */
    cStrand *strand = new cStrand(Apto::String("abcdxyzabcdxyzabcd"));
    EXPECT_TRUE(strand);
    EXPECT_EQ(strand->AsString(), Apto::String("abcdxyzabcdxyzabcd"));
    EXPECT_EQ(strand->GetSize(), sizeof("abcdxyzabcdxyzabcd")-1);
    delete strand;
  }
};
