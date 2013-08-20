/*
 *  cStrand.cc
 *  Avida
 *
 *  File created by Kaben Nanlohy on 7/31/13.
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
 */

#include "cStrand.h"
#include "apto/core/Algorithms.h"
#include "apto/scheduler/Util.h"

#include <iostream>
#include <cmath>

using namespace std;


Apto::Scheduler::ProbabilisticDynamic::~ProbabilisticDynamic() { }

void Apto::Scheduler::ProbabilisticDynamic::AdjustPriority(int entry_id, double priority) {
  m_index.SetWeight(entry_id, priority);
}

int Apto::Scheduler::ProbabilisticDynamic::Next() {
  if (m_index.TotalWeight() == 0) { return -1; }
  return m_index.FindPosition(m_rng->GetDouble(m_index.TotalWeight()));
}


Apto::Scheduler::ProbabilisticDynamic::DynamicWeightedIndex::DynamicWeightedIndex(int size)
: m_size(size), m_item_weight(size), m_subtree_weight(size)
{
  m_item_weight.SetAll(0.0);
  m_subtree_weight.SetAll(0.0);
}


void Apto::Scheduler::ProbabilisticDynamic::DynamicWeightedIndex::SetWeight(int entry_id, double weight) {
  m_item_weight[entry_id] = weight;
  while (true) {
    const int left_id = leftChildOf(entry_id);
    const int right_id = rightChildOf(entry_id);
    const double left_subtree = (left_id >= m_size) ? 0.0 : m_subtree_weight[left_id];
    const double right_subtree = (right_id >= m_size) ? 0.0 : m_subtree_weight[right_id];
    m_subtree_weight[entry_id] = m_item_weight[entry_id] + left_subtree + right_subtree;
    if (entry_id == 0) { break; }
    entry_id = parentOf(entry_id);
  }
}

int Apto::Scheduler::ProbabilisticDynamic::DynamicWeightedIndex::findPosition(double position, int root_id) {
  assert(position < m_subtree_weight[root_id]);
  // First, see if we should just return this node.
  if (position < m_item_weight[root_id]) { return root_id; }
  // If not, then see if we should search in the left subtree...
  position -= m_item_weight[root_id];
  const int left_id = leftChildOf(root_id);
  assert (left_id < m_size);
  if (position < m_subtree_weight[left_id]) { return findPosition(position, left_id); }
  // Otherwise we must look in the right subtree...
  position -= m_subtree_weight[left_id];
  const int right_id = rightChildOf(root_id);
  assert (right_id < m_size);
  assert (position < m_subtree_weight[right_id]);
  return findPosition(position, right_id);
}

void Apto::Scheduler::ProbabilisticDynamic::DynamicWeightedIndex::Resize(int size) {
  if (size < m_size) { for (int i = m_size - 1; i >= size; i--) { SetWeight(i, 0.); } }
  m_item_weight.Resize(size, 0.);
  m_subtree_weight.Resize(size, 0.);
  m_size = size;
}


void Apto::Scheduler::ProbabilisticID::AdjustIDPriority(int id, double priority) {
  int pos = -1;
  if (!id2pos.Has(id)) {
    pos = id2pos.GetSize();
    if (GetSize() <= pos) {
      Resize(pos + 1);
      pos2id.Push(id);
    } else { pos2id[pos] = id; }
    id2pos[id] = pos;
  } else {
    pos = id2pos[id];
  }
  AdjustPriority(pos, priority);
}

int Apto::Scheduler::ProbabilisticID::NextID() { return pos2id[Next()]; }


Apto::Scheduler::IntegratedDynamic::IntegratedDynamic(int entry_count) : m_priority_chart(entry_count) {
  m_priority_chart.SetAll(0.0);
}

Apto::Scheduler::IntegratedDynamic::~IntegratedDynamic() {
  for (int i = 0; i < m_node_array.GetSize(); i++) { if (m_node_array[i] != NULL) { delete m_node_array[i]; } }
}

void Apto::Scheduler::IntegratedDynamic::AdjustPriority(int entry_id, double priority) {
  double old_priority = m_priority_chart[entry_id];
  if (old_priority == priority) { return; }
  m_priority_chart[entry_id] = priority;
  // Re-adjust the lists.
  Util::Priority old_p_comp(old_priority);
  Util::Priority new_p_comp(priority);
  int priority_magnitude = Max(old_p_comp.NumBits(), new_p_comp.NumBits());
  for (int i = 0; i < priority_magnitude; i++) {
    bool old_bit = old_p_comp.Bit(i);
    bool new_bit = new_p_comp.Bit(i);
    if (old_bit && !new_bit) {
      // Remove the item from this node...
      m_node_array[i]->Remove(entry_id);
      if (m_node_array[i]->size == 0) { removeNode(i); }
    }
    if (!old_bit && new_bit) {
      // Add the item from this node...
      if (i >= m_node_array.GetSize() || !m_node_array[i]) { insertNode(i); }
      m_node_array[i]->Insert(entry_id);
    }
  }
}

int Apto::Scheduler::IntegratedDynamic::Next() {
  assert(m_node_array.GetSize() > 0);  // Running scheduler w/ no entries!
  const int last_id = m_node_array.GetSize() - 1;
  // Make sure there are entries in the scheduler!
  if (m_node_array[last_id] == NULL) { return -1; }
  int next_id = -1;
  while (next_id < 0) { next_id = m_node_array[last_id]->Next(); }
  return next_id;
}

void Apto::Scheduler::IntegratedDynamic::Resize(int size) {
  if (size < GetSize()) for (int i = GetSize() - 1; i >= size; i--) { AdjustPriority(i, 0.); }
  for (int i = 0; i < m_node_array.GetSize(); i++) { m_node_array[i]->active_array.Resize(size, 0); }
  m_priority_chart.Resize(size, 0.);
}


void Apto::Scheduler::IntegratedDynamic::insertNode(int node_id) {
  // Test if trying to create node that already exists.
  assert(node_id >= m_node_array.GetSize() || m_node_array[node_id] == NULL);
  Node* new_node = new Node(GetSize(), node_id);
  if (node_id >= m_node_array.GetSize()) { resizeNodes(node_id); }
  m_node_array[node_id] = new_node;
  // Find the node to mark as the 'prev'.
  for (int prev_id = node_id + 1; prev_id < m_node_array.GetSize(); prev_id++) {
    Node* prev_node = m_node_array[prev_id];
    if (prev_node) {
      new_node->prev = prev_node;
      prev_node->next = new_node;
      prev_node->process_size = (1 << (prev_id - node_id - 1));
      break;
    }
  }
  // And find the node to mark as the 'next'.
  for (int next_id = node_id - 1; next_id >= 0; next_id--) {
    Node* next_node = m_node_array[next_id];
    if (next_node) {
      new_node->next = next_node;
      next_node->prev = new_node;
      new_node->process_size = (1 << (node_id - next_id - 1));
      break;
    }
  }
}

void Apto::Scheduler::IntegratedDynamic::removeNode(int node_id) {
  assert(m_node_array[node_id] != NULL); // Trying to remove non-existant node.
  Node* old_node = m_node_array[node_id];
  Node* next_node = old_node->next;
  Node* prev_node = old_node->prev;
  m_node_array[node_id] = NULL;
  if (next_node) { next_node->prev = prev_node; }
  if (prev_node) {
    prev_node->next = next_node;
    prev_node->process_size = old_node->process_size * prev_node->process_size * 2;
  }
  if (node_id == m_node_array.GetSize() - 1) {
    if (!old_node->next) { resizeNodes(0); }
    else { resizeNodes(old_node->next->node_id); }
  }
  delete old_node;
}

void Apto::Scheduler::IntegratedDynamic::resizeNodes(int new_max)
{
  int old_size = m_node_array.GetSize();
  int new_size = new_max + 1;  // 0 to new_max...
  // Clean up tail portions of the array being cut off.
  for (int i = new_size; i < old_size; i++) { if (m_node_array[i]) { delete m_node_array[i]; } }
  m_node_array.Resize(new_size);
  // Mark as NULL any new cells added to the array.
  for (int i = old_size; i < new_size; i++) { m_node_array[i] = NULL; }
}

void Apto::Scheduler::IntegratedDynamic::Node::Insert(int item_id) {
  assert(item_id >= 0 && item_id < active_array.GetSize());  // Illegal ID
  // If this item is already active in this node, ignore this call...
  if (active_array[item_id] != 0) { return; }
  // See if we're dealing with a new first_entry...
  if (first_entry == -1 || item_id < first_entry) {
    active_array[item_id] = first_entry;
    first_entry = item_id;
  }
  else {
    // Otherwise find the predecessor to this item in the list...
    int prev_item;
    for (prev_item = item_id - 1; prev_item >= 0; prev_item--) { if (active_array[prev_item] != 0) break; }
    assert(prev_item >= 0);  // prev_item is first, but not identified.
    // Make the predecessor point to it, and have it point to the CPU that
    // the old predecessor pointed to.
    active_array[item_id] = active_array[prev_item];
    active_array[prev_item] = item_id;
  }
  size++;
}

void Apto::Scheduler::IntegratedDynamic::Node::Remove(int item_id)
{
  assert(item_id >= 0 && item_id < active_array.GetSize()); // Illegal ID
  // If this item is already inactive, ignore this call...
  if (active_array[item_id] == 0) { return; }
  // If this is the first_entry, adjust it!
  if (first_entry == item_id) { first_entry = active_array[item_id]; }
  else {
    // Find the predecessor to this item in the list...
    int prev_item;
    for (prev_item = item_id - 1; prev_item >= 0; prev_item--) { if (active_array[prev_item] != 0) break; }
    assert(prev_item >= 0);  // prev_item is first, but not identified.
    // Make the predecessor point to the item removed used to point to.
    active_array[prev_item] = active_array[item_id];
  }
  active_array[item_id] = 0;
  size--;
}

// Execute everything on list, and then shift to calling the next node.
// Wait for the next node to return a -1 before shifting back to this one.
int Apto::Scheduler::IntegratedDynamic::Node::Next() {
  // Alternate between this node's Process and the next's.
  if (execute == false) {
    // If there is a next node, we may be working on it...
    int next_id = -1;
    if (next != NULL) { next_id = next->Next(); }
    // If next_id is a -1, either we don't have a next node, or else it
    // is finished with its execution.
    if (next_id == -1) {
      execute = true;
      process_count = 0;
      active_entry = -1;
    }
    return next_id;
  }
  // Find the next active_entry...
  // If we were at the end of the list, start over...
  if (active_entry == -1) active_entry = first_entry;
  // If this entry no longer exists, hunt for the next active entry manually...
  else if (active_array[active_entry] == 0) {
    while (active_entry < active_array.GetSize() && active_array[active_entry] == 0) { active_entry++; }
    if (active_entry == active_array.GetSize()) { active_entry = -1; }
  }
  // Otherwise, if the entry does exist, we can just look the next one up.
  else active_entry = active_array[active_entry];
  // If we have now hit the end of this list, move on to the next node.
  if (active_entry == -1) {
    process_count++;
    if (process_count >= process_size) { execute = false; }
  }
  return active_entry;
}


void Apto::Scheduler::IntegratedID::AdjustIDPriority(int id, double priority) {
  int pos = -1;
  if (!id2pos.Has(id)) {
    pos = id2pos.GetSize();
    if (GetSize() <= pos) {
      Resize(pos + 1);
      pos2id.Push(id);
    } else { pos2id[pos] = id; }
    id2pos[id] = pos;
  } else {
    pos = id2pos[id];
  }
  AdjustPriority(pos, priority);
}

int Apto::Scheduler::IntegratedID::NextID() { return pos2id[Next()]; }




int IDMgr::NextID() {
  int id = -1;
  if (0 < m_recycled_ids.GetSize()) {
    /* Try to recycle old IDs from things that used to be in the index but were removed. */
    id = m_recycled_ids.Pop();
  } else {
    /* Otherwise assign a new ID. */
    id = m_max_ct;
    m_max_ct += 1;
  }
  m_obj_ct += 1;
  return id;
}
void IDMgr::RecycleID(int id) {
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
int cBasicLabelUtils::SeqLen(const int id) const {
  const int base = GetNumNops();
  int seqlen = 0, offset = 0, next_offset = 1;
  while (next_offset <= id) {
    offset = next_offset;
    next_offset = (offset * base) + 1;
    seqlen++;
  }
  return seqlen;
}
int cBasicLabelUtils::Rotate(const int id, const int rot) const {
  /* First convert to a rotated sequence. */
  const int base = GetNumNops();
  int seqlen = 0, offset = 0, next_offset = 1;
  while (next_offset <= id) {
    offset = next_offset;
    next_offset = (offset * base) + 1;
    seqlen++;
  }
  int seqval = id - offset;
  int seq[seqlen];
  for (int i = seqlen-1; 0 <= i; i--) {
    /* Rotate happens here. */
    seq[i] = (seqval + rot) % base;
    seqval /= base;
  }
  /* Now convert rotated sequence to a new ID. */
  int new_id = 0;
  for (int i = 0; i < seqlen; i++) {
    new_id *= base;
    new_id += seq[i] + 1;
  }
  return new_id;
}
int cBasicLabelUtils::ReverseComplement(const int id, const int rot) const {
  /* First convert to a reverse-complement sequence. */
  const int base = GetNumNops();
  int seqlen = 0, offset = 0, next_offset = 1;
  while (next_offset <= id) {
    offset = next_offset;
    next_offset = (offset * base) + 1;
    seqlen++;
  }
  int seqval = id - offset;
  int seq[seqlen];
  for (int i = seqlen-1; 0 <= i; i--) {
    /* Reverse complement happens here. */
    seq[seqlen - 1 - i] = (seqval + rot) % base;
    seqval /= base;
  }
  /* Now convert reverse-complement sequence to a new ID. */
  int new_id = 0;
  for (int i = 0; i < seqlen; i++) {
    new_id *= base;
    new_id += seq[i] + 1;
  }
  return new_id;
}


bool cLabelIdx::Has(int id) { return m_id2idx.Has(id); }
cLabel* cLabelIdx::Get(int id) { return (Has(id))?(m_objs[m_id2idx[id]]):((cLabel*)(0)); }
cLabel* cLabelIdx::Insert(int id) {
  if (Has(id)) { return m_objs[m_id2idx[id]]; }
  int idx = NextID();
  if (m_objs.GetSize() < GetMaxCt()) { m_objs.Resize(GetMaxCt()); }
  cLabel* ptr = new cLabel;
  ptr->m_id = id;
  m_objs[idx] = ptr;
  m_id2idx[id] = idx;
  return ptr;
}
bool cLabelIdx::Delete(int id) {
  if (Has(id)) {
    delete m_objs[m_id2idx[id]];
    m_objs[m_id2idx[id]] = (cLabel*)(0);
    RecycleID(m_id2idx[id]);
    m_id2idx.Remove(id);
    return true;
  }
  return false;
}
cLabelIdx::Iterator::Iterator(cLabelIdx& idx) : m_idx(idx), m_it(m_idx.m_id2idx.Keys()) { ; }
int cLabelIdx::Iterator::ID() {
  if (m_it.Get()) { return *m_it.Get(); }
  return -1;
}
cLabel* cLabelIdx::Iterator::Get() {
  if (m_it.Get()) { return m_idx.m_objs[m_idx.m_id2idx[*m_it.Get()]]; }
  return (cLabel*)(0);
}
cLabel* cLabelIdx::Iterator::Next() {
  if (m_it.Next()) { return m_idx.m_objs[m_idx.m_id2idx[*m_it.Get()]]; }
  return (cLabel*)(0);
}


bool cSeqIdx::Has(int id) { return (0 <= id) && (id < GetMaxCt()) && (m_objs[id] != (cSequence*)(0)); }
bool cSeqIdx::Has(const Apto::String &str) { return m_str2id.Has(str); }
cSequence* cSeqIdx::Get(int id) { return (Has(id))?(m_objs[id]):((cSequence*)(0)); }
int cSeqIdx::GetID(const Apto::String &str) { return (Has(str))?(m_str2id.Get(str)):(-1); }
Apto::String cSeqIdx::GetString(int id) { return (Has(id))?(m_id2str.Get(id)):(""); }
cSequence* cSeqIdx::Insert(const Apto::String &str) {
  if (Has(str)) { return m_objs[GetID(str)]; }
  int id = NextID();
  if (m_objs.GetSize() < GetMaxCt()) { m_objs.Resize(GetMaxCt()); }
  cSequence* ptr = new cSequence;
  ptr->m_id = id;
  m_objs[id] = ptr;
  m_str2id[str] = id;
  m_id2str[id] = str;
  return ptr;
}
bool cSeqIdx::Delete(int id) {
  if (Has(id)) {
    m_str2id.Remove(m_id2str[id]);
    m_id2str.Remove(id);
    delete m_objs[id];
    m_objs[id] = (cSequence*)(0);
    RecycleID(id);
    return true;
  }
  return false;
}
bool cSeqIdx::Delete(const Apto::String &str) { return Delete(GetID(str)); }


int cFSMDB::InsertSequence(const Apto::String& seq) {
  if (!m_seqs.Has(seq)) {
    Apto::Array<cHit, Apto::Smart> hits(bScanForLabels(seq, m_label_utils));
    cSequence* seq_ptr = m_seqs.Insert(seq);
    int seq_id = seq_ptr->ID();
    Apto::Map<int, Apto::Set<int> > lbl_sites;
    for (int i=0; i<hits.GetSize(); i++) {
      const int lbl_id = hits[i].Lbl();
      const int lbl_pos = hits[i].Pos();
      cLabel* lbl_ptr(m_lbls.Insert(lbl_id));
      lbl_ptr->m_seq_ids.Insert(seq_id);
      lbl_sites[lbl_id].Insert(lbl_pos);
    }
    for (Apto::Map<int, Apto::Set<int> >::KeyIterator it = lbl_sites.Keys(); it.Next();) {
      int lbl_id = *it.Get();
      seq_ptr->m_lbl_sites[lbl_id].Resize(lbl_sites[lbl_id].GetSize());
      int j=0;
      for (Apto::Set<int>::Iterator pos_it = lbl_sites[lbl_id].Begin(); pos_it.Next();) {
        seq_ptr->m_lbl_sites[lbl_id][j] = *pos_it.Get();
        j++;
      }
    }
    return seq_id;
  } else { return m_seqs.GetID(seq); }
}
void cFSMDB::UnlinkSeqLbls(int seq_id) {
  cSequence* seq_ptr = m_seqs.Get(seq_id);
  if (seq_ptr) {
    for (Apto::Map<int, Apto::Array<int> >::KeyIterator it = seq_ptr->m_lbl_sites.Keys(); it.Next();) {
      int lbl_id = *it.Get();
      cLabel* lbl_ptr(m_lbls.Get(lbl_id));
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
  /* Get or create strand and sequence */
  cStrand* strand_ptr = m_bindables.Create<cStrand>();
  int strand_id = strand_ptr->ID();
  int seq_id = InsertSequence(seq);
  cSequence* seq_ptr = m_seqs.Get(seq_id);
  /* Link strand and sequence. */
  strand_ptr->m_seq_id = seq_id;
  seq_ptr->m_strand_ids.Insert(strand_id);
  /*
  Schedule collisions with priority proportional to strand length.

  This is a place where we might be able to control the probability that two
  particular molecules will collide.
  */
  m_collision_scheduler.Resize(m_bindables.GetMaxCt());
  m_collision_scheduler.AdjustPriority(strand_id, seq.GetSize());
  return strand_id;
}
bool cFSMDB::RemoveStrand(int strand_id) {
  cStrand* strand_ptr = m_bindables.Get<cStrand>(strand_id);
  if (strand_ptr) {
    cSequence* seq_ptr = m_seqs.Get(strand_ptr->m_seq_id);
    seq_ptr->m_strand_ids.Remove(strand_id);
    if (seq_ptr->m_strand_ids.GetSize() < 1) { RemoveSequence(strand_ptr->m_seq_id); }
    /* Unschedule collisions for this strand ID. */
    m_collision_scheduler.AdjustPriority(strand_id, 0.);
    return m_bindables.Delete(strand_id);
  } else { return false; }
}


Apto::Array<cHit, Apto::Smart> bScanForLabels(const Apto::String& seq, const cBasicLabelUtils& label_utils) {
  /*
   Scan along the strand looking for label sections (consisting of certain
   Nops). Whenever we're in a label section, make a note of the position of
   every possible label.
   */
  Apto::Array<cHit, Apto::Smart> hits;
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
        hits.Push(cHit(id, pos));
      }
    }
  }
  return hits;
}

