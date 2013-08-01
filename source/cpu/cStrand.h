/*
 *  cStrand.h
 *  Avida
 *
 *  Copyright 2013 Michigan State University. All rights reserved.
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
 */

#ifndef cStrand_h
#define cStrand_h

#include "cCPUMemory.h"
#include "tList.h"

class cStrand;
class cStrandLabelPositions;

/*
 A strand is at heart just a string. But we're going to add per-strand label
 tracking. A label is a substring consisting of certain special characters.
 */
class cStrand : public cCPUMemory {
public:
  Apto::Map<int, Apto::Array<int> > m_label_ids_by_position;
  Apto::Map<int, cStrandLabelPositions > m_label_positions_by_id;
public:
  cStrand(const cStrand& in_strand) : cCPUMemory(in_strand) {}
  cStrand(const cCPUMemory& in_memory) : cCPUMemory(in_memory) {}
  cStrand(const InstructionSequence& in_genome) : cCPUMemory(in_genome) {}
  cStrand(const Apto::String& in_string) : cCPUMemory(in_string) {}
  explicit cStrand(int size = 1) : cCPUMemory(size) {}
  ~cStrand() { ; }
};

/*
 There can be several copies of each label in a given strand, each copy
 appearing at a different position. For each label ID, the corresponding
 cStrandLabelPositions records the locations of each copy of the label.
 
 The strand's label index is a hashmap of these cStrandLabelPositions records,
 and is keyed by the distinct label IDs.
 */
class cStrandLabelPositions {
public:
  Apto::Array<int> m_label_positions;
  /*
   m_strand_index_node is initially set to NULL. When the owning strand is
   placed in the main index, backrefs into the index are stored in
   m_strand_index_node, one backref for each label-ID. The backrefs are
   recorded so we can quickly find the parts of the index that need updating
   when we want to remove the strand from the index.
   */
  tListNode<cStrand> *m_index_backref;
public:
  cStrandLabelPositions(int num_positions = 0)
  : m_label_positions(num_positions)
  , m_index_backref(0)
  {}
};

#endif
