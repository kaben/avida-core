/*
 *  cStrand.cc
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

#include "cStrand.h"

#include <iostream>

using namespace std;


int IDBase::NextID() {
  int id = -1;
  if (0 < m_recycled_ids.GetSize()) {
    /* Try to recycle old IDs from things that used to be in the index but were removed. */
    id = m_recycled_ids.Pop();
  } else {
    /* Otherwise assign a new ID. */
    id = m_next_new_id;
    m_next_new_id += 1;
  }
  m_obj_ct += 1;
  return id;
}
void IDBase::RecycleID(int id) {
  m_recycled_ids.Push(id);
  m_obj_ct -= 1;
}


int cBasicLabelUtils::IsNop(char c) const {
  c -= m_zero_char;
  return (0 <= c) && (c < m_num_nops);
}
int cBasicLabelUtils::Seq2ID(const Apto::String& seq) const {
  const int base = GetNumNops();
  const int zero_char = m_zero_char;
  const int seqlen = seq.GetSize();
  int id = 0;
  for (int i = 0; i < seqlen; i++) {
    id *= base;
    id += seq[i] - zero_char + 1;
  }
  return id;
}
Apto::StringBuffer cBasicLabelUtils::ID2Seq(const int id) const {
  const int base = GetNumNops();
  const int zero_char = m_zero_char;
  int seqlen = 0, offset = 0, next_offset = 1;
  while (next_offset <= id) {
    offset = next_offset;
    next_offset = (offset * base) + 1;
    seqlen++;
  }
  int seqval = id - offset;
  Apto::StringBuffer seq(seqlen);
  for (int i = seqlen-1; 0 <= i; i--) {
    seq[i] = zero_char + seqval % base;
    seqval /= base;
  }
  return seq;
}


bool cLabelIdx::Has(int id) { return m_id2obj.Has(id); }
cLabel* cLabelIdx::Get(int id) { return (Has(id))?(&m_id2obj.Get(id)):(0); }
int cLabelIdx::GetSize() { return m_id2obj.GetSize(); }
cLabel* cLabelIdx::Insert(int id) {
  cLabel* ptr = 0;
  if (Has(id)) {
    ptr = Get(id);
  } else {
    ptr = &m_id2obj.Get(id);
    ptr->m_id = id;
  }
  return ptr;
}
bool cLabelIdx::Delete(int id) {
  if (Has(id)) {
    m_id2obj.Remove(id);
    return true;
  }
  return false;
}
cLabelIdx::IDIter cLabelIdx::IDs() { return m_id2obj.Keys(); }
cLabelIdx::ObjIter cLabelIdx::Objs() { return m_id2obj.Values(); }


bool cSeqIdx::Has(int id) { return m_id2obj.Has(id); }
bool cSeqIdx::Has(const Apto::String &str) { return m_str2id.Has(str); }
cSequence* cSeqIdx::Get(int id) { return (Has(id))?(&m_id2obj.Get(id)):(0); }
int cSeqIdx::GetID(const Apto::String &str) { return (Has(str))?(m_str2id.Get(str)):(-1); }
Apto::String cSeqIdx::GetString(int id) { return (Has(id))?(m_id2str.Get(id)):(""); }
cSequence* cSeqIdx::Insert(const Apto::String &str) {
  int id = GetID(str);
  if (id < 0) {
    id = NextID();
    m_str2id[str] = id;
    m_id2str[id] = str;
    m_id2obj.Get(id).m_id = id;
  }
  return &m_id2obj.Get(id);
}
bool cSeqIdx::Delete(int id) {
  if (Has(id)) {
    m_str2id.Remove(m_id2str[id]);
    m_id2str.Remove(id);
    m_id2obj.Remove(id);
    RecycleID(id);
    return true;
  }
  return false;
}
bool cSeqIdx::Delete(const Apto::String &str) { return Delete(GetID(str)); }
cSeqIdx::IDIter cSeqIdx::IDs() { return m_id2obj.Keys(); }
cSeqIdx::ObjIter cSeqIdx::Objs() { return m_id2obj.Values(); }


void cFSMDB::LinkSeqLblPos(int seq_id, int lbl_id, int lbl_pos) {
  cSequence* seq_ptr = m_seqs.Get(seq_id);
  if (seq_ptr) {
    cLabel* lbl_ptr = m_lbls.Insert(lbl_id);
    if (lbl_ptr) {
      lbl_ptr->m_seq_ids.Insert(seq_id);
      seq_ptr->m_lbl_sites[lbl_id].Insert(lbl_pos);
    }
  }
}
int cFSMDB::InsertSequence(const Apto::String& seq) {
  if (!m_seqs.Has(seq)) {
    Apto::Array<cLabelHit, Apto::Smart> hits(bScanForLabels(seq, m_label_utils));
    int seq_id = m_seqs.Insert(seq)->ID();
    for (int i=0; i<hits.GetSize(); i++) { LinkSeqLblPos(seq_id, hits[i].Lbl(), hits[i].Pos()); }
    return seq_id;
  } else { return m_seqs.GetID(seq); }
}
void cFSMDB::UnlinkSeqLbls(int seq_id) {
  cSequence* seq_ptr = m_seqs.Get(seq_id);
  if (seq_ptr) {
    for (Apto::Map<int, Apto::Set<int> >::KeyIterator it = seq_ptr->m_lbl_sites.Keys(); it.Next();) {
      int lbl_id = *it.Get();
      cLabel* lbl_ptr = m_lbls.Get(lbl_id);
      lbl_ptr->m_seq_ids.Remove(seq_id);
      if (lbl_ptr->m_seq_ids.GetSize() < 1) { m_lbls.Delete(lbl_id); }
    }
    seq_ptr->m_lbl_sites.Clear();
  }
}
bool cFSMDB::RemoveSequence(int seq_id) {
  UnlinkSeqLbls(seq_id);
  return m_seqs.Delete(seq_id);
}
int cFSMDB::CreateStrand(const Apto::String &seq) {
  cStrand* strand_ptr = m_strands.Create();
  int seq_id = InsertSequence(seq);
  strand_ptr->m_seq_id = seq_id;
  int strand_id = strand_ptr->ID();
  m_seqs.Get(seq_id)->m_strand_ids.Insert(strand_id);
  return strand_id;
}
bool cFSMDB::RemoveStrand(int strand_id) {
  cStrand* strand_ptr = m_strands.Get(strand_id);
  if (strand_ptr) {
    cSequence* seq_ptr = m_seqs.Get(strand_ptr->m_seq_id);
    seq_ptr->m_strand_ids.Remove(strand_id);
    if (seq_ptr->m_strand_ids.GetSize() < 1) { RemoveSequence(strand_ptr->m_seq_id); }
    return m_strands.Delete(strand_id);
  } else { return false; }
}


Apto::Array<cLabelHit, Apto::Smart> bScanForLabels(const Apto::String& seq, const cBasicLabelUtils& label_utils) {
  /*
   Scan along the strand looking for label sections (consisting of certain
   Nops). Whenever we're in a label section, make a note of the position of
   every possible label.
   */
  Apto::Array<cLabelHit, Apto::Smart> hits;
  const int default_max_size = label_utils.GetMaxLabelSize();
  int label_section_start = 0, label_section_size = 0;
  bool in_label_section = false;
  for (int pos = 0; pos < seq.GetSize(); pos++) {
    if (in_label_section){
      /* Check whether we're exiting label section. */
      if (!label_utils.IsNop(seq[pos])) in_label_section = false;
      else label_section_size++;
    } else {
      /* Check whether we're entering label section. */
      if (label_utils.IsNop(seq[pos])){
        in_label_section = true;
        label_section_start = pos;
        label_section_size = 1;
      }
    }
    if (in_label_section){
      int max_size = default_max_size;
      if (label_section_size < max_size) max_size = label_section_size;
      for (int label_size = 1; label_size <= max_size; label_size++) {
        int pos = label_section_start + label_section_size - label_size;
        int id = label_utils.Seq2ID(seq.Substring(pos, label_size));
        hits.Push(cLabelHit(id, pos));
      }
    }
  }
  return hits;
}

