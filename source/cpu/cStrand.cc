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

#include <algorithm>
#include <iostream>
#include <cmath>
#include <vector>

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
  const int left_id = leftChildOf(root_id); assert (left_id < m_size);
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
  for (int i = 0; i < m_node_array.GetSize(); i++) {
    if (m_node_array[i] != NULL) {
      m_node_array[i]->active_array.Resize(size, 0);
    }
  }
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

void cHalfBinding::Print() {
  cout << "cHalfBinding::Print(): ID: " << ID() << ", m_parent_id: " << m_parent_id << ", m_lbl_id: " << m_lbl_id << ", m_lbl_pos: " << m_lbl_pos << ", m_lbl_len: " << m_lbl_len << ", m_other_half_binding_id: " << m_other_half_binding_id << endl;
}

Apto::Array<int> cBindable::GetBindings(cFSMDB &db) {
  Apto::Set<int> half_binding_set;
  for (Apto::Map<int, Apto::Set<int> >::KeyIterator kit = m_bindpts.Keys(); kit.Next();) {
    for (Apto::Set<int>::Iterator it = m_bindpts[*kit.Get()].Begin(); it.Next();) {
       half_binding_set.Insert(*it.Get());
    }
  }
  Apto::Array<int> half_binding_ary;
  for (Apto::Set<int>::Iterator it = half_binding_set.Begin(); it.Next();) {
    half_binding_ary.Push(*it.Get());
  }
  return half_binding_ary;
}

Apto::Map<int, Apto::Array<int> > &cStrand::GetLabels(cFSMDB &db) {
  return db.m_seqs.Get(m_seq_id)->m_labels;
}
Apto::String cStrand::AsString(cFSMDB &db) {
  return db.m_seqs.GetString(m_seq_id);
}


Apto::Map<int, Apto::Array<int> > &cFSM::GetLabels(cFSMDB &db) {
  return db.m_fsm_defs.Get(m_fsm_def_id)->m_labels;
}
Apto::String cFSM::AsString(cFSMDB &db) {
  return "";
}

int cNFA::Transition(int symbol_id, cFSMDB &db) {
  cNFADef *nfa_def = db.m_fsm_defs.Get<cNFADef>(m_fsm_def_id);
  cout << "symbol_id: " << symbol_id << ", current m_state_id: " << m_state_id;
  if (nfa_def != NULL) {
    if (nfa_def->m_transition_relation.Has(m_state_id) && nfa_def->m_transition_relation[m_state_id].Has(symbol_id)) {
      int num_states = nfa_def->m_transition_relation[m_state_id][symbol_id].GetSize();
      if (1 < num_states) {
        cout << ", num_states: " << num_states;
        int state_idx = m_rng->GetInt(0, num_states);
        cout << ", state_idx: " << state_idx;
        m_state_id = nfa_def->m_transition_relation[m_state_id][symbol_id][state_idx];
      } else {
        m_state_id = nfa_def->m_transition_relation[m_state_id][symbol_id][0];
      }
    }
    if (nfa_def->m_function_relation.Has(m_state_id)) {
      int num_functors = nfa_def->m_function_relation[m_state_id].GetSize();
      cout << ", num_functors: " << num_functors;
      int functor_id = -1;
      if (1 < num_functors) {
        int functor_idx = m_rng->GetInt(0, num_functors);
        cout << ", functor_idx: " << functor_idx;
        functor_id = nfa_def->m_function_relation[m_state_id][functor_idx];
        cout << ", functor_id: " << functor_id;
      } else {
        functor_id = nfa_def->m_function_relation[m_state_id][0];
      }
      //FSMFunctor functor;
      cout << endl;
      db.m_functors[functor_id](ID());
    }
  }
  cout << ", next m_state_id: " << m_state_id << endl;
  return m_state_id;
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


cFSMFunctorObject::cFSMFunctorObject(cFSMDB &db)
: m_db(db)
{}
void cFSMFunctorObject::Function0(int caller_id){
  cout << "cFSMFunctorObject::Function0(" << caller_id << ")" << endl;
}
void cFSMFunctorObject::Function1(int caller_id){
  cout << "cFSMFunctorObject::Function1(" << caller_id << ")" << endl;
}
void cFSMFunctorObject::Function2(int caller_id){
  cout << "cFSMFunctorObject::Function2(" << caller_id << ")" << endl;
}
void cFSMFunctorObject::Function3(int caller_id){
  cout << "cFSMFunctorObject::Function3(" << caller_id << ")" << endl;
}
void cFSMFunctorObject::Function4(int caller_id){
  cout << "cFSMFunctorObject::Function4(" << caller_id << ")" << endl;
}
void cFSMFunctorObject::Function5(int caller_id){
  cout << "cFSMFunctorObject::Function5(" << caller_id << ")" << endl;
}
void cFSMFunctorObject::Function6(int caller_id){
  cout << "cFSMFunctorObject::Function6(" << caller_id << ")" << endl;
}
void cFSMFunctorObject::Function7(int caller_id){
  cout << "cFSMFunctorObject::Function7(" << caller_id << ")" << endl;
}
void cFSMFunctorObject::Function8(int caller_id){
  cout << "cFSMFunctorObject::Function8(" << caller_id << ")" << endl;
}
void cFSMFunctorObject::Function9(int caller_id){
  cout << "cFSMFunctorObject::Function9(" << caller_id << ")" << endl;
}


cFSMDB::cFSMDB(int rng_seed)
: m_rng(new Apto::RNG::AvidaRNG(rng_seed))
//, m_kinetics(m_rng)
, m_collision_scheduler(0, m_rng)
, m_unbinding_scheduler(0)
, m_functors(10)
, m_functor_object(*this)
{
  m_functors[0] = FSMFunctor(&m_functor_object, &cFSMFunctorObject::Function0);
  m_functors[1] = FSMFunctor(&m_functor_object, &cFSMFunctorObject::Function1);
  m_functors[2] = FSMFunctor(&m_functor_object, &cFSMFunctorObject::Function2);
  m_functors[3] = FSMFunctor(&m_functor_object, &cFSMFunctorObject::Function3);
  m_functors[4] = FSMFunctor(&m_functor_object, &cFSMFunctorObject::Function4);
  m_functors[5] = FSMFunctor(&m_functor_object, &cFSMFunctorObject::Function5);
  m_functors[6] = FSMFunctor(&m_functor_object, &cFSMFunctorObject::Function6);
  m_functors[7] = FSMFunctor(&m_functor_object, &cFSMFunctorObject::Function7);
  m_functors[8] = FSMFunctor(&m_functor_object, &cFSMFunctorObject::Function8);
  m_functors[9] = FSMFunctor(&m_functor_object, &cFSMFunctorObject::Function9);
}
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
      seq_ptr->m_labels[lbl_id].Resize(lbl_sites[lbl_id].GetSize());
      int j=0;
      for (Apto::Set<int>::Iterator pos_it = lbl_sites[lbl_id].Begin(); pos_it.Next();) {
        seq_ptr->m_labels[lbl_id][j] = *pos_it.Get();
        j++;
      }
    }
    return seq_id;
  } else { return m_seqs.GetID(seq); }
}
void cFSMDB::UnlinkSeqLbls(int seq_id) {
  cSequence* seq_ptr = m_seqs.Get(seq_id);
  if (seq_ptr) {
    for (Apto::Map<int, Apto::Array<int> >::KeyIterator it = seq_ptr->m_labels.Keys(); it.Next();) {
      int lbl_id = *it.Get();
      cLabel* lbl_ptr(m_lbls.Get(lbl_id));
      lbl_ptr->m_seq_ids.Remove(seq_id);
      if (lbl_ptr->m_seq_ids.GetSize() < 1) { m_lbls.Delete(lbl_id); }
    }
    seq_ptr->m_labels.Clear();
  }
}
bool cFSMDB::RemoveSequence(int seq_id) {
  UnlinkSeqLbls(seq_id);
  return m_seqs.Delete(seq_id);
}
int cFSMDB::CreateStrand() {
  /* Get or create strand and sequence */
  cStrand* strand_ptr = m_bindables.Create<cStrand>();
  strand_ptr->m_seq_id = -1;
  int strand_id = strand_ptr->ID();
  m_collision_scheduler.Resize(m_bindables.GetMaxCt());
  m_collision_scheduler.AdjustPriority(strand_id, 0);
  return strand_id;
}
int cFSMDB::CreateStrand(const Apto::String &seq) {
  int strand_id = CreateStrand();
  AssociateSeqToStrand(strand_id, seq);
  return strand_id;
}
void cFSMDB::AssociateSeqToStrand(int strand_id, const Apto::String &seq) {
  cStrand* strand_ptr = m_bindables.Get<cStrand>(strand_id); assert(NULL != strand_ptr);
  int old_seq_id = strand_ptr->m_seq_id;

  /* Disassociate old sequence, if any. */
  cSequence* old_seq_ptr(NULL);
  if (0 <= old_seq_id) {
    old_seq_ptr = m_seqs.Get(old_seq_id); assert(NULL != old_seq_ptr);
    /* Unlink strand and old sequence. */
    old_seq_ptr->m_strand_ids.Remove(strand_id);
    strand_ptr->m_seq_id = -1;
  }

  /* Associate new sequence, if any. */
  if (0 < seq.GetSize()) {
    int new_seq_id = InsertSequence(seq);
    cSequence* new_seq_ptr = m_seqs.Get(new_seq_id); assert(NULL != new_seq_ptr);
    /* Link strand and new sequence. */
    strand_ptr->m_seq_id = new_seq_id;
    new_seq_ptr->m_strand_ids.Insert(strand_id);
  }

  /* If old sequence no longer has any associated strands, delete it. */
  if (0 <= old_seq_id) {
    if (old_seq_ptr->m_strand_ids.GetSize() < 1) { RemoveSequence(old_seq_id); }
  }
  m_collision_scheduler.AdjustPriority(strand_id, seq.GetSize());
}
void cFSMDB::RemoveStrand(int strand_id) {
  /* Extract a list of half binding IDs in the parent.  */
  cStrand* par = m_bindables.Get<cStrand>(strand_id); assert(NULL != par);
  Apto::Set<int> par_hb_ids(CollapseSetMap(par->m_bindpts));
  /* Delete all half bindings; also find and delete their other halves.  */
  for (Apto::Set<int>::Iterator it = par_hb_ids.Begin(); it.Next();) {
    int hb_id = *it.Get();
    cHalfBinding *hb = m_half_bindings.Get(hb_id);
    int pos = hb->m_lbl_pos;
    int len = hb->m_lbl_len;

    int o_hb_id = hb->m_other_half_binding_id;
    cHalfBinding *o_hb = m_half_bindings.Get(o_hb_id);
    int o_pos = o_hb->m_lbl_pos;
    cBindable *o_par = m_bindables.Get(o_hb->m_parent_id);
    for (int i=0; i< len; i++) {
      o_par->m_bindpts[o_pos + i].Remove(o_hb_id);
      if (0 == o_par->m_bindpts[o_pos+i].GetSize()) { o_par->m_bindpts.Remove(o_pos+i); }
    }
    m_half_bindings.Delete(o_hb_id);
    m_half_bindings.Delete(hb_id);
  }
  par->m_bindpts.Clear();
  /* Disassociate strand and sequence. */
  AssociateSeqToStrand(strand_id, "");
  /* Delete strand. */
  m_bindables.Delete(strand_id);
}
void cFSMDB::SplitStrand(int strand_id, int at_pos, int &ret_d0_id, int &ret_d1_id) {
  cStrand* par = m_bindables.Get<cStrand>(strand_id); assert(NULL != par);
  Apto::String par_str(par->AsString(*this));
  int par_len = par_str.GetSize(); assert(at_pos <= par_len);
  ret_d0_id = CreateStrand(par_str.Substring(0, at_pos));
  ret_d1_id = CreateStrand(par_str.Substring(at_pos));
  cStrand *d0 = m_bindables.Get<cStrand>(ret_d0_id);
  cStrand *d1 = m_bindables.Get<cStrand>(ret_d1_id);

  /*
  Convert daughters' pos-(label-array) maps to pos-(label-set) maps. The latter
  will allow us to use Has() to check whether a label appears at the given
  position in the daughter.
  */
  Apto::Map<int, Apto::Set<int> > d0_lbls(AsSetMap(d0->GetLabels(*this)));
  Apto::Map<int, Apto::Set<int> > d1_lbls(AsSetMap(d1->GetLabels(*this)));
  /* Extract a list of half binding IDs from the parent. */
  Apto::Set<int> par_hb_ids(CollapseSetMap(par->m_bindpts));
  /* Remove half bindings from parent. */
  par->m_bindpts.Clear();
  /* Transfer half bindings to daughters. */
  for (Apto::Set<int>::Iterator it = par_hb_ids.Begin(); it.Next();) {
    int hb_id = *it.Get();
    cHalfBinding *hb = m_half_bindings.Get(hb_id);
    int id = hb->m_lbl_id;
    int pos = hb->m_lbl_pos;
    int len = hb->m_lbl_len;
    bool found = false;
    if (pos < at_pos) {
      /* Transfer to daughter 0. */
      if (d0_lbls.Has(id) && d0_lbls[id].Has(pos)) {
        found = true;
        hb->m_parent_id = ret_d0_id;
        for (int i=0; i < len; i++) { d0->m_bindpts[pos+i].Insert(hb_id); }
      }
    } else {
      /* Transfer to daughter 1. */
      int ofs_pos = pos - at_pos;
      if (d1_lbls.Has(id) && d1_lbls[id].Has(ofs_pos)) {
        found = true;
        hb->m_parent_id = ret_d1_id;
        hb->m_lbl_pos = ofs_pos;
        for (int i=0; i < len; i++) { d1->m_bindpts[ofs_pos+i].Insert(hb_id); }
      }
    }
    /* Delete bindings that couldn't be transferred. */
    if (!found) {
      int o_hb_id = hb->m_other_half_binding_id;
      cHalfBinding *o_hb = m_half_bindings.Get(o_hb_id);
      int o_pos = o_hb->m_lbl_pos;
      cBindable *o_par = m_bindables.Get(o_hb->m_parent_id);
      for (int i=0; i< len; i++) {
        o_par->m_bindpts[o_pos+i].Remove(o_hb_id);
        if (0 == o_par->m_bindpts[o_pos+i].GetSize()) { o_par->m_bindpts.Remove(o_pos+i); }
      }
      m_half_bindings.Delete(o_hb_id);
      m_half_bindings.Delete(hb_id);
    }
  }
  /* Delete parent strand. */
  RemoveStrand(strand_id);
}
void cFSMDB::JoinStrands(int strand_0_id, int strand_1_id, int &ret_daughter_id) {
  cStrand* p0 = m_bindables.Get<cStrand>(strand_0_id); assert(NULL != p0);
  cStrand* p1 = m_bindables.Get<cStrand>(strand_1_id); assert(NULL != p1);
  Apto::String p0_str(p0->AsString(*this)), p1_str(p1->AsString(*this));
  int p0_len = p0_str.GetSize(), p1_len = p1_str.GetSize();
  ret_daughter_id = CreateStrand(p0_str+p1_str);
  cStrand *d = m_bindables.Get<cStrand>(ret_daughter_id);
  /*
  Convert daughter's pos-(label-array) map to pos-(label-set) map. The latter
  will allow us to use Has() to check whether a label appears at the given
  position in the daughter.
  */
  Apto::Map<int, Apto::Set<int> > d_lbls(AsSetMap(d->GetLabels(*this)));
  /* Extract lists of half binding IDs from the parents. */
  Apto::Set<int> p0_hb_ids(CollapseSetMap(p0->m_bindpts));
  Apto::Set<int> p1_hb_ids(CollapseSetMap(p1->m_bindpts));
  /* Remove half bindings from parents. */
  p0->m_bindpts.Clear(); p1->m_bindpts.Clear();
  /* Transfer first parent's half bindings to daughter. */
  for (Apto::Set<int>::Iterator it = p0_hb_ids.Begin(); it.Next();) {
    int hb_id = *it.Get();
    cHalfBinding *hb = m_half_bindings.Get(hb_id);
    int id = hb->m_lbl_id;
    int pos = hb->m_lbl_pos;
    int len = hb->m_lbl_len;
    hb->m_parent_id = ret_daughter_id;
    for (int i=0; i < len; i++) { d->m_bindpts[pos+i].Insert(hb_id); }
  }
  /* Transfer second parent's half bindings to daughter. */
  for (Apto::Set<int>::Iterator it = p1_hb_ids.Begin(); it.Next();) {
    int hb_id = *it.Get();
    cHalfBinding *hb = m_half_bindings.Get(hb_id);
    int id = hb->m_lbl_id;
    int pos = hb->m_lbl_pos + p0_len; /* Offset by first parent's length. */
    int len = hb->m_lbl_len;
    hb->m_parent_id = ret_daughter_id;
    for (int i=0; i < len; i++) { d->m_bindpts[pos+i].Insert(hb_id); }
  }
  /* Delete parent strands. */
  RemoveStrand(strand_0_id);
  RemoveStrand(strand_1_id);
}
int cFSMDB::CreateFSMBootstrap() {
  cFSMBootstrap* ptr = m_bindables.Create<cFSMBootstrap>();
  int id = ptr->ID();
  return 0;
}
bool cFSMDB::SingleCollision() {
  /* Get two molecules to collide. */
  int bindable_id_0 = m_collision_scheduler.Next();
  int bindable_id_1 = m_collision_scheduler.Next();
  if ((bindable_id_0 == -1) || (bindable_id_1 == -1)) { return false; }
  return Collide(bindable_id_0, bindable_id_1);
}
bool cFSMDB::Collide(int bbl_id_0, int bbl_id_1) {
  cBindable* bbl_0 = m_bindables.Get<cBindable>(bbl_id_0);
  cBindable* bbl_1 = m_bindables.Get<cBindable>(bbl_id_1);
  /* Brainstorming ... display bindables as strings. */
  Apto::String str_0(bbl_0->AsString(*this));
  Apto::String str_1(bbl_1->AsString(*this));
  cout << "bbl_id_0:" << bbl_id_0 << ":str_0:" << str_0 << ", bbl_id_1:" << bbl_id_1 << ":str_1:" << str_1 << endl;
  Apto::Map<int, Apto::Array<int> > &labels_0(bbl_0->GetLabels(*this)), &labels_1(bbl_1->GetLabels(*this));
  /*
  (Here we're just choosing the shorter of the two sequences in hopes that
  we can reduce processing time a bit.)
  */
  int bbl_id_a = -1, bbl_id_b = -1;
  cBindable *bbl_a, *bbl_b;
  Apto::Map<int, Apto::Array<int> > *labels_a(0), *labels_b(0);
  int lbl_ct_0 = labels_0.GetSize();
  int lbl_ct_1 = labels_1.GetSize();
  if (lbl_ct_0 < lbl_ct_1) {
    bbl_a = bbl_0; bbl_b = bbl_1;
    bbl_id_a = bbl_id_0; bbl_id_b = bbl_id_1;
    labels_a = &labels_0; labels_b = &labels_1;
  } else {
    bbl_a = bbl_1; bbl_b = bbl_0;
    bbl_id_a = bbl_id_1; bbl_id_b = bbl_id_0;
    labels_a = &labels_1; labels_b = &labels_0;
  }
  /* Walk through possible bindings; organize by length. */
  Apto::Array<Apto::Array<int, Apto::Smart> > label_ids_by_length(m_label_utils.GetMaxLabelSize() + 1);
  for (Apto::Map<int, Apto::Array<int> >::KeyIterator it = labels_a->Keys(); it.Next();) {
    const int lbl_id = *it.Get();
    const int lbl_len = m_label_utils.SeqLen(lbl_id);
    const int rvc_id = m_label_utils.ReverseComplement(lbl_id);
    if (labels_b->Has(rvc_id)) { label_ids_by_length[lbl_len].Push(lbl_id); }
  }
  /* Start trying to bind, giving preferences to longer binding sites. */
  for (int len = label_ids_by_length.GetSize() - 1; 0 < len; len--) {
    /*
    Catalog all possible bindings, searching first by labels of this
    length, together with their complements; and then by positions of each
    label and its complement.
    */
    std::vector<int> halfbinding_ids;
    int lbl_ct = label_ids_by_length[len].GetSize();
    for (int j = 0; j < lbl_ct; j++) {
      /*
      Gather info about each possible pairing of a label and its
      complement.
      */
      int lbl_id = label_ids_by_length[len][j];
      int rvc_id = m_label_utils.ReverseComplement(lbl_id);
      int lbl_position_ct = (*labels_a)[lbl_id].GetSize();
      int rvc_position_ct = (*labels_b)[rvc_id].GetSize();
      for (int k = 0; k < lbl_position_ct; k++) {
        for (int l = 0; l < rvc_position_ct; l++) {
          cHalfBinding *lbl_hb = m_half_bindings.Create();
          cHalfBinding *rvc_hb = m_half_bindings.Create();
          m_unbinding_scheduler.Resize(m_half_bindings.GetMaxCt());
          lbl_hb->Set(bbl_id_a, lbl_id, (*labels_a)[lbl_id][k], len, rvc_hb->ID());
          rvc_hb->Set(bbl_id_b, rvc_id, (*labels_b)[rvc_id][l], len, lbl_hb->ID());
          halfbinding_ids.push_back(lbl_hb->ID());
        }
      }
    }
    /* Shuffle binding ids. */
    /* TODO: add vector shuffling to Avida, based on Avida's RNG. */
    std::random_shuffle(halfbinding_ids.begin(), halfbinding_ids.end());
    /*
    Randomly select bindings, and discard bindings that won't work because
    their positions have already been used, or bindings that randomly fail
    to bind...
    */
    /*
    TODO: establish probability distributions for bindings.

    I chose these probabilities pretty much arbitrarily. For the purposes
    of brainstorming, I made longer bindings more probable than shorter
    bindings. But this isn't necessarily how I think binding will
    eventually work.

    Update: switched point binding probability to 0.10910128185966073, so
    probability of 6-point binding is about 0.5 (or
    1-(1-0.10910128185966073)**6). (Computed as 1-math.exp(math.log(0.5)/6.).)
    */
    double point_binding_probability = 0.10910128185966073;
    double binding_probability = 1. - pow(1. - point_binding_probability, len);
    for (std::vector<int>::iterator it=halfbinding_ids.begin(); it!=halfbinding_ids.end(); ++it) {
      /* Get the two halves of the binding candidate. */
      int half_binding_id = *it;
      cHalfBinding *lbl_hb = m_half_bindings.Get(half_binding_id);
      int other_half_binding_id = lbl_hb->m_other_half_binding_id;
      cHalfBinding *rvc_hb = m_half_bindings.Get(other_half_binding_id);
      /*
      Check to see whether the pair can bind -- that is, whether their
      respective binding positions are available for binding.
      */
      bool can_bind = true;
      for (int j = 0; j < len; j++) {
        int lbl_chk = lbl_hb->m_lbl_pos + j;
        int rvc_chk = rvc_hb->m_lbl_pos + j;
        /*
        For brainstorming, we're only allowing one binding at a particular
        point, but code for the class allows more than one binding.

        I think that whether to allow more than one binding depends upon the
        reactants; but that most of the time only one binding would be allowed.
        An exception could be needed if, for example, we want certain molecules
        to be able to catalyze an unbinding reaction, in which case the
        catalyst would bind to an already-bound position, and then the double
        binding would immediately completely unbind.

        Giving a precise quantitative description of this kind of reaction will
        be easier once we have chemical kinetics.
        */
        if (bbl_a->m_bindpts.Has(lbl_chk) && (bbl_a->m_bindpts[lbl_chk].GetSize()>0)){ can_bind = false; }
        if (bbl_b->m_bindpts.Has(rvc_chk) && (bbl_b->m_bindpts[rvc_chk].GetSize()>0)){ can_bind = false; }
      }
      /* Flip a loaded coin to see whether the binding would succeed. */
      bool does_bind = (m_rng->GetDouble(0., 1.) < binding_probability);
      if (can_bind && does_bind) {
        /* If binding is possible and will succeed, bind! */
        /* Brainstorming ... display binding info. */
        cout << "  binding: " << m_label_utils.ID2Seq(lbl_hb->m_lbl_id) << " (" << lbl_hb->m_lbl_pos << "), " << m_label_utils.ID2Seq(rvc_hb->m_lbl_id) << " (" << rvc_hb->m_lbl_pos << ")" << endl;
        for (int j = 0; j < len; j++) {
          /* In each molecule, mark all bound positions. */
          int lbl_chk = lbl_hb->m_lbl_pos + j;
          int rvc_chk = rvc_hb->m_lbl_pos + j;
          bbl_a->m_bindpts[lbl_chk].Insert(lbl_hb->ID());
          bbl_b->m_bindpts[rvc_chk].Insert(rvc_hb->ID());
          /* Brainstorming ... display binding info. */
          cout << "    " << lbl_hb->m_lbl_pos + j << ":" << rvc_hb->m_lbl_pos + j << endl;
        }
        m_unbinding_scheduler.AdjustPriority(lbl_hb->ID(), len);
        m_unbinding_scheduler.AdjustPriority(rvc_hb->ID(), len);
      } else {
        /* Delete failed bindings from database. */
        m_unbinding_scheduler.AdjustPriority(lbl_hb->ID(), 0.);
        m_unbinding_scheduler.AdjustPriority(rvc_hb->ID(), 0.);
        m_half_bindings.Delete(lbl_hb->ID());
        m_half_bindings.Delete(rvc_hb->ID());
      }
    }
  }
  /* FIXME: currently always returns true; should give meaning to retval, or remove. */
  return true;
}
bool cFSMDB::SingleUnbinding() {
  int lbl_hb_id = m_unbinding_scheduler.Next();
  if (lbl_hb_id == -1) { return false; }
  return Unbind(lbl_hb_id);
}
bool cFSMDB::Unbind(int lbl_hb_id) {
  cHalfBinding *lbl_hb = m_half_bindings.Get(lbl_hb_id);
  int rvc_hb_id = lbl_hb->m_other_half_binding_id;
  cHalfBinding *rvc_hb = m_half_bindings.Get(rvc_hb_id);
  /*
  TODO: establish probability distributions for bindings.

  I chose these probabilities pretty much arbitrarily. For the purposes of
  brainstorming, I made longer bindings more probable than shorter bindings.
  But this isn't necessarily how I think binding will eventually work.

  Update: switched point binding probability to 0.10910128185966073, so
  probability of 6-point binding is about 0.5 (or
  1-(1-0.10910128185966073)**6).
  */
  double point_binding_probability = 0.10910128185966073;
  double unbinding_probability = pow(1. - point_binding_probability, lbl_hb->m_lbl_len);
  bool does_unbind = (m_rng->GetDouble(0., 1.) < unbinding_probability);
  if (does_unbind) {
    cout << "unbinding: " << m_label_utils.ID2Seq(lbl_hb->m_lbl_id) << " (" << lbl_hb->m_lbl_pos << "), " << m_label_utils.ID2Seq(rvc_hb->m_lbl_id) << " (" << rvc_hb->m_lbl_pos << ")" << endl;
    int bbl_id_a = lbl_hb->m_parent_id;
    int bbl_id_b = rvc_hb->m_parent_id;
    cBindable* bbl_a = m_bindables.Get<cBindable>(bbl_id_a);
    cBindable* bbl_b = m_bindables.Get<cBindable>(bbl_id_b);
    Apto::String str_a(bbl_a->AsString(*this));
    Apto::String str_b(bbl_b->AsString(*this));
    cout << "  bbl_id_a:" << bbl_id_a << ":str_a:" << str_a << ", bbl_id_b:" << bbl_id_b << ":str_b:" << str_b << endl;
    for (int j = 0; j < lbl_hb->m_lbl_len; j++) {
      /* In each molecule, mark all bound positions. */
      int lbl_chk = lbl_hb->m_lbl_pos + j;
      int rvc_chk = rvc_hb->m_lbl_pos + j;
      bbl_a->m_bindpts[lbl_chk].Remove(lbl_hb->ID());
      bbl_b->m_bindpts[rvc_chk].Remove(rvc_hb->ID());
      /* Brainstorming ... display binding info. */
      cout << "  " << lbl_hb->m_lbl_pos + j << ":" << rvc_hb->m_lbl_pos + j << endl;
    }
    m_unbinding_scheduler.AdjustPriority(lbl_hb->ID(), 0.);
    m_unbinding_scheduler.AdjustPriority(rvc_hb->ID(), 0.);
    m_half_bindings.Delete(lbl_hb->ID());
    m_half_bindings.Delete(rvc_hb->ID());
  }
  /* FIXME: currently always returns true; should give meaning to retval, or remove. */
  return true;
}
bool cFSMDB::SingleRebinding() {
  int lbl_hb_id = m_unbinding_scheduler.Next();
  if (lbl_hb_id == -1) { return false; }
  cHalfBinding *lbl_hb = m_half_bindings.Get(lbl_hb_id);
  int rvc_hb_id = lbl_hb->m_other_half_binding_id;
  cHalfBinding *rvc_hb = m_half_bindings.Get(rvc_hb_id);
  int bbl_id_a = lbl_hb->m_parent_id;
  int bbl_id_b = rvc_hb->m_parent_id;
  cout << "rebinding: " << m_label_utils.ID2Seq(lbl_hb->m_lbl_id) << " (" << lbl_hb->m_lbl_pos << "), " << m_label_utils.ID2Seq(rvc_hb->m_lbl_id) << " (" << rvc_hb->m_lbl_pos << ")" << endl;
  Unbind(lbl_hb_id);
  Collide(bbl_id_a, bbl_id_b);
  return true;
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

