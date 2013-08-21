/*
 *  unittests/core/Strand.cc
 *  avida-core
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
 *  Authors: Kaben G. Nanlohy <kaben.nanlohy@gmail.com>
 *
 */

#include "cpu/cStrand.h"

#include "apto/scheduler/Util.h"

#include "gtest/gtest.h"

#include <iostream>
#include <cstdlib>      // std::rand, std::srand

using namespace std;


/*
BOOKMARK 20130819-2217
- Possible solution to need for overlapping bindings implemented by changing
  cBindable::m_bindpts from Map<int, int> to Map<int, Array<int> >. Thus each
  bindable position has not a possibly-empty slot for a single half-binding id,
  but a possibly-empty slot for a stack of half-binding ids.
BOOKMARK 20130819-1303
Thoughts:
- Possible processing steps:
  - ProcessChemicalKinetics()
  - ProcessFSMs()
- It seems to me this should just be a ProcessStep() function, which executes
  any of the following according to some schedule:
  - chemical kinetics
  - FSMs
- It occurs to me that FSMs might be thought of as in the domain of chemical
  kinetics.
  - The problem with this is that I'm thinking of kinetics in terms of
    collisions on the one hand, and in terms of FSM steps on the other hand. To
    put them in the same terms would mean adding possibly-avoidable complexity.
BOOKMARK 20130819-0800
Thoughts:
- When two strands are bound, other strands can't bind at the bound positions.
  - This is in part intentional, as it can be used for regulation of gene
    expression, which is its intended purpose.
  - On the other hand, it seems like certain FSMs should be able to temporarily
    unbind and then rebind the bound strands, by analogy to RNA synthesis from
    DNA.
  - So a bound pair of strands needs to be able to collide with an FSM, which
    should be able to unbind the pair.
  - Perhaps multiple bindings should be allowed, although a second binding
    might have a lower probability.
  - I might also permit a binding to disrupt another.
  - This might all be abstracted into a chemical kinetics manager.
Next steps:
- Resume work on collisions.
*/


//class cFSMdb;
//
//class cStrandLibrary;
//class cStrand;
//
///* Using Apto::String instead. */
////class bString;
//class cSequence;
//
//class bBindingSiteLibrary;
//class bBindingSite;
//class bBindingSiteID;
//class bBindingSitePosition;
//
//class cFSMLibrary;
//class cFSM;
//class cFSMProgramBase;
//class cFSMProgramStrandCPU;
//class cFSMProgramStrandBootstrap;
//class cFSMProgramStrandEncoding;
//
//class cHit;
//class bBindingLibrary;
//class bBinding;
//class bBindingCriteria;
//
//Apto::Array<cHit, Apto::Smart> bScanForLabels(int seq_id, cFSMdb &db);
//
//class cSequence {
//public:
//  Apto::String m_str;
//  Apto::Map<int, Apto::Set<int> > m_lbl_sites;
//  typedef typename Apto::Map<int, Apto::Set<int> >::KeyIterator LabelIter;
//  typedef typename Apto::Set<int>::Iterator PositionIter;
//public:
//  cSequence()
//  {}
//public:
//  void InsertLabelPos(int lbl_id, int lbl_pos) { m_lbl_sites[lbl_id].Insert(lbl_pos); }
//  bool GetLabelHitCt(int lbl_id) { return m_lbl_sites[lbl_id].GetSize(); }
//  bool HasLabel(int lbl_id) { return m_lbl_sites.Has(lbl_id); }
//  bool RemoveLabelPos(int lbl_id, int lbl_pos) { return m_lbl_sites[lbl_id].Remove(lbl_pos); }
//  bool RemoveLabel(int lbl_id) { return m_lbl_sites.Remove(lbl_id); }
//  Apto::String GetString() { return m_str; }
//  LabelIter Labels() { return m_lbl_sites.Keys(); }
//  PositionIter Positions(int label_id) { return m_lbl_sites[label_id].Begin(); }
//};
//
//class cLabel {
//protected:
//  Apto::Set<int> m_seq_ids;
//  Apto::Set<int> m_fsm_def_ids;
//public:
//  cLabel()
//  {}
//public:
//  void InsertSeq(int seq_id) { m_seq_ids.Insert(seq_id); }
//  bool HasSeq(int seq_id) { return m_seq_ids.Has(seq_id); }
//  int GetSeqCount() { return m_seq_ids.GetSize(); }
//  bool RemoveSeq(int seq_id) { m_seq_ids.Remove(seq_id); }
//
//  void InsertFSMDef(int fsm_def_id) { m_fsm_def_ids.Insert(fsm_def_id); }
//  bool HasFSMDef(int fsm_def_id) { return m_fsm_def_ids.Has(fsm_def_id); }
//  bool GetFSMDefCount() { return m_fsm_def_ids.GetSize(); }
//  bool RemoveFSMDef(int fsm_def_id) { m_fsm_def_ids.Remove(fsm_def_id); }
//};
//
//class cLabelDeletemeIdx {
//protected:
//  Apto::Map<int, cLabel> m_id2obj;
//public:
//  typedef typename Apto::Map<int, cLabel> AM;
//  typedef typename AM::KeyIterator IDIter;
//  typedef typename AM::ValueIterator ObjIter;
//public:
//  bool Has(int id) { return m_id2obj.Has(id); }
//  cLabel* Get(int id) { return (Has(id))?(&m_id2obj.Get(id)):(0); }
//  int GetSize() { return m_id2obj.GetSize(); }
//  cLabel* Insert(int id) {
//    cLabel* ptr = 0;
//    if (Has(id)) {
//      ptr = Get(id);
//    } else {
//      ptr = &m_id2obj.Get(id);
//    }
//    return ptr;
//  }
//  bool Delete(int id) {
//    if (Has(id)) {
//      m_id2obj.Remove(id);
//      return true;
//    }
//    return false;
//  }
//  IDIter IDs() { return m_id2obj.Keys(); }
//  ObjIter Objs() { return m_id2obj.Values(); }
//};
//
//class cStrand {
//public:
//  int m_seq_id;
//public:
//  cStrand()
//  : m_seq_id(-1)
//  {}
//public:
//};
//
//
//
///* bSeqIDIdx stores unique strings, and associates an ID with each. */
//class bSeqIDIdx : public IDBase {
//protected:
//  Apto::Map<int, cSequence> m_id2obj;
//  Apto::Map<Apto::String, int> m_str2id;
//  const Apto::String m_nil;
//public:
//  /* Constructor. Optional arg nil is returned to indicate id not found in index when GetStr(id) is called. */
//  bSeqIDIdx(const Apto::String nil = ""): m_nil(nil) {}
//public:
//  bool Has(int id) { return m_id2obj.Has(id); }
//  bool Has(const Apto::String &str) { return m_str2id.Has(str); }
//  Apto::SmartPtr<cSequence> Get(int id) { return (Has(id))?(&m_id2obj.Get(id)):(0); }
//  int GetID(const Apto::String &str) { return (Has(str))?(m_str2id.Get(str)):(-1); }
//  Apto::String GetString(int id) {
//    /* Indicate if id isn't in the index by returning m_no_such_str. */
//    Apto::SmartPtr<cSequence> ptr = Get(id);
//    if (ptr) { return ptr->GetString(); }
//    return m_nil;
//  }
//  /* Add str to index if it isn't there. Return corresponding ID. */
//  int Insert(const Apto::String &str) {
//    int id = GetID(str);
//    if (id < 0) {
//      id = NextID();
//      Apto::SmartPtr<cSequence> ptr = &m_id2obj.Get(id);
//      m_str2id[str] = id;
//
//      ptr->m_str = str;
//    }
//    return id;
//  }
//  /* If str is in index, remove and return corresponding ID. Otherwise return -1. */
//  int Delete(const Apto::String &str) {
//    int id = GetID(str);
//    if (0 <= id) {
//      /* If 0 <= id then str is in the index, so remove it. */
//      m_id2obj.Remove(id);
//      m_str2id.Remove(str);
//      RecycleID(id);
//    }
//    return id;
//  }
//  bool Delete(int id) {
//    Apto::SmartPtr<cSequence> ptr = Get(id);
//    if (ptr) {
//      m_str2id.Remove(ptr->m_str);
//      m_id2obj.Remove(id);
//      RecycleID(id);
//      return true;
//    }
//    return false;
//  }
//};
//
//
//class bSeqLabelAssoc {
//public:
//  int m_label_id;
//  int m_seq_id;
//  int m_seq_pos;
//public:
//  bSeqLabelAssoc()
//  : m_label_id(-1)
//  , m_seq_id(-1)
//  , m_seq_pos(-1)
//  {}
//public:
//};
//
//
//class bSeqLabelIdx {
//public:
//  IntAryMap m_labels_by_seq;
//  IntAryMap m_seqs_by_label;
//public:
//public:
//};


/*
cFSMdb owns all objects, and each object has an ID.

*/



//class cFSMdb {
//public:
//  bSeqIDIdx m_seq_idx;
//  ObjDeleteMeIdx<cStrand> m_strand_idx;
//  cLabelDeletemeIdx m_label_idx;
//
//  cBasicLabelUtils m_label_utils;
//public:
//  cStrand* GetStrand(int id);
//
//  int InsertSequence(const Apto::String& sequence); 
//  Apto::String GetSequenceString(int id);
//  bool HasSequence(int id);
//  bool HasSequence(Apto::SmartPtr<cSequence> &ptr);
//  bool HasSequence(const Apto::String &str);
//  cSequence *GetSequence(int id); 
//  int GetSequenceID(cSequence *ptr); 
//  int GetSequenceID(const Apto::String &str); 
//  bool RemoveSequence(const Apto::String& sequence);
//  bool RemoveSequence(Apto::SmartPtr<cSequence> &ptr);
//  bool RemoveSequence(int id);
//
//  cLabel* InsertLabel(int id); 
//  bool HasLabel(cLabel* ptr);
//  bool HasLabel(int id);
//  cLabel* GetLabel(int id);
//  int GetLabelID(cLabel* ptr); 
//  int LabelRefCt(cLabel* ptr);
//  int LabelRefCt(int id);
//public:
//  void LinkSeqLblPos(int seq_id, int lbl_id, int lbl_pos);
//  void UnlinkSeqLblPos(int seq_id, int lbl_id, int lbl_pos);
//  void UnlinkSeqLbl(int seq_id, int lbl_id);
//  void UnlinkSeqLbls(int seq_id);
//protected:
//  bool RemoveLabel(cLabel* ptr);
//  bool RemoveLabel(int id);
//};
//
//cStrand* cFSMdb::GetStrand(int id) { return m_strand_idx.Get(id); }
//
//int cFSMdb::InsertSequence(const Apto::String& seq) {
//  int seq_id = -1;
//  if (m_seq_idx.Has(seq)) {
//    /* Already known sequence; don't scan known sequences. */
//    seq_id = m_seq_idx.Insert(seq);
//  } else {
//    /* Scan new sequence for labels. */
//    seq_id = m_seq_idx.Insert(seq);
//    Apto::Array<cHit, Apto::Smart> hits(bScanForLabels(seq_id, *this));
//    for (int i=0; i<hits.GetSize(); i++) {
//      LinkSeqLblPos(seq_id, hits[i].Lbl(), hits[i].Pos());
//    }
//  }
//  return seq_id;
//}
//
//bool cFSMdb::HasSequence(int id) { return m_seq_idx.Has(id); }
//bool cFSMdb::HasSequence(Apto::SmartPtr<cSequence> ptr) { return m_seq_idx.Has(ptr); }
//bool cFSMdb::HasSequence(const Apto::String &str) { return m_seq_idx.Has(str); }
//
//void cFSMdb::UnlinkSeqLblPos(int seq_id, int lbl_id, int lbl_pos) {
//  Apto::SmartPtr<cSequence> seq_ptr = GetSequence(seq_id);
//  if (seq_ptr) {
//    seq_ptr->RemoveLabelPos(lbl_id, lbl_pos);
//    if (seq_ptr->GetLabelHitCt(lbl_id) <= 0) {
//      seq_ptr->RemoveLabel(lbl_id);
//      cLabel* lbl_ptr = GetLabel(lbl_id);
//      if (lbl_ptr) {
//        lbl_ptr->RemoveSeq(seq_id);
//        if (LabelRefCt(lbl_ptr) < 1) { RemoveLabel(lbl_ptr); }
//      }
//    }
//  }
//}
//void cFSMdb::UnlinkSeqLbl(int seq_id, int lbl_id) {
//  Apto::SmartPtr<cSequence> seq_ptr = GetSequence(seq_id);
//  cLabel* lbl_ptr = GetLabel(lbl_id);
//  if (seq_ptr && lbl_ptr) {
//    lbl_ptr->RemoveSeq(seq_id);
//    seq_ptr->RemoveLabel(lbl_id);
//    if (LabelRefCt(lbl_ptr) < 1) { RemoveLabel(lbl_ptr); }
//  }
//}
//void cFSMdb::UnlinkSeqLbls(int seq_id) {
//  Apto::SmartPtr<cSequence> seq_ptr = GetSequence(seq_id);
//  if (seq_ptr) {
//    for (cSequence::LabelIter it = seq_ptr->Labels(); it.Next();) {
//      int lbl_id = *it.Get();
//      cLabel* lbl_ptr = GetLabel(lbl_id);
//      lbl_ptr->RemoveSeq(seq_id);
//      /* Do not do this, as the Map iterators don't handle deletion during iteration:
//      seq_ptr->RemoveLabel(lbl_id);
//      */
//      if (LabelRefCt(lbl_ptr) < 1) { RemoveLabel(lbl_ptr); }
//    }
//    /* Do this instead of deleting while iterating. */
//    seq_ptr->m_lbl_sites.Clear();
//  }
//}
//
//Apto::String cFSMdb::GetSequenceString(int id) { return m_seq_idx.GetString(id); }
//Apto::SmartPtr<cSequence> cFSMdb::GetSequence(int id) { return m_seq_idx.Get(id); }
//int cFSMdb::GetSequenceID(cSequence *ptr) { return m_seq_idx.GetID(ptr); }
//int cFSMdb::GetSequenceID(const Apto::String &str) { return m_seq_idx.GetID(str); }
//bool cFSMdb::RemoveSequence(const Apto::String& str) {
//  RemoveSequence(GetSequenceID(str));
//}
//bool cFSMdb::RemoveSequence(Apto::SmartPtr<cSequence> ptr) {
//  RemoveSequence(GetSequenceID(ptr));
//}
//bool cFSMdb::RemoveSequence(int id) {
//  UnlinkSeqLbls(id);
//  m_seq_idx.Delete(id);
//}
//
//cLabel* cFSMdb::InsertLabel(int id) { return m_label_idx.Insert(id); }
//bool cFSMdb::HasLabel(cLabel* ptr) { return m_label_idx.Has(ptr); }
//bool cFSMdb::HasLabel(int id) { return m_label_idx.Has(id); }
//cLabel* cFSMdb::GetLabel(int id) { return m_label_idx.Get(id); }
//int cFSMdb::GetLabelID(cLabel* ptr) { return m_label_idx.GetID(ptr); }
//int cFSMdb::LabelRefCt(cLabel* ptr) {
//  if (m_label_idx.Has(ptr)) { return ptr->GetSeqCount(); }
//  return -1;
//}
//int cFSMdb::LabelRefCt(int id) {
//  if (m_label_idx.Has(id)) { return m_label_idx.Get(id)->GetSeqCount(); }
//  return -1;
//}
//bool cFSMdb::RemoveLabel(cLabel* ptr) { return m_label_idx.Delete(ptr); }
//bool cFSMdb::RemoveLabel(int id) { return m_label_idx.Delete(id); }


//class bBindingSitePosition {
//public:
//};
//
//class cFSM {
//public:
//};
//
//class cFSMProgramBase {
//public:
//};
//
//class cFSMProgramStrandCPU {
//public:
//};
//
//class cFSMProgramStrandBootstrap {
//public:
//};
//
//class cFSMProgramStrandEncoding {
//public:
//};
//
//class bBinding {
//public:
//};
//
//class bBindingCriteria{
//public:
//};



TEST(cSeqIdx, brainstorm){
  cSeqIdx idx;

  EXPECT_EQ(0, idx.Insert(Apto::String("blah"))->ID());
  EXPECT_EQ(0, idx.GetID("blah"));
  EXPECT_EQ("blah", idx.GetString(0));

  /* Handles nonsense input? */
  EXPECT_EQ(-1, idx.GetID("(nonsense)"));
  EXPECT_EQ("", idx.GetString(-1));

  /* Try to insert same thing a second time; id shouldn't change. */
  EXPECT_EQ(0, idx.Insert(Apto::String("blah"))->ID());
  EXPECT_EQ(0, idx.GetID("blah"));

  EXPECT_EQ(1, idx.Insert(Apto::String("foo"))->ID());
  EXPECT_EQ(1, idx.GetID("foo"));
  EXPECT_EQ("foo", idx.GetString(1));

  EXPECT_TRUE(idx.Delete(Apto::String("blah")));
  /* Try to lookup something that was removed; should indicate nonsense. */
  EXPECT_EQ(-1, idx.GetID("blah"));

  /* Try to remove something that has already been removed; should indicate nonsene. */
  EXPECT_FALSE(idx.Delete(Apto::String("blah")));

  /* id 0 should now be available for reuse. */
  EXPECT_EQ(0, idx.Insert(Apto::String("bar"))->ID());
  EXPECT_EQ(0, idx.GetID("bar"));
  EXPECT_EQ("bar", idx.GetString(0));

  /* There should be no more ids available for reuse, so the next should be 2. */
  EXPECT_EQ(2, idx.Insert(Apto::String("fubar"))->ID());
  EXPECT_EQ(2, idx.GetID("fubar"));
  EXPECT_EQ("fubar", idx.GetString(2));

  EXPECT_TRUE(idx.Delete(Apto::String("fubar")));
  EXPECT_EQ(-1, idx.GetID("fubar"));
  EXPECT_FALSE(idx.Delete(2));
}

namespace cLabelIdxTests {

  TEST(cLabelIdx, instantiation){
    cLabelIdx idx;
    /* Index should initially be empty, so have zero size! */
    EXPECT_EQ(0, idx.GetSize());
  }

  TEST(cLabelIdx, handle_get_nonsense){
    cLabelIdx idx;
    /* Try to get using a nonsense ID; should return null ptr. */
    EXPECT_FALSE(idx.Get(-1));
    /* Try to get using a bad ID; should return null ptr. */
    EXPECT_FALSE(idx.Get(0));
    /* Get with bad ID should not accidentally add to index! */
    EXPECT_FALSE(idx.Get(0));
    EXPECT_EQ(0, idx.GetSize());
  }

  TEST(cLabelIdx, insert){
    cLabelIdx idx;
    /* Create should return non-null ptr. */
    cLabel* ptr = idx.Insert(5);
    EXPECT_TRUE(idx.Insert(5));
    /* Object should have ID 5! */
    EXPECT_EQ(5, ptr->ID());
    EXPECT_EQ(1, idx.GetSize());
  }

  TEST(cLabelIdx, persistence){
    cLabelIdx idx;
    /* Changes to a managed object should persist. */
    int id = -1;
    {
      cLabel* ptr_0 = idx.Insert(5);
      EXPECT_EQ(0, ptr_0->m_seq_ids.GetSize());
      EXPECT_FALSE(ptr_0->m_seq_ids.Has(7));
      ptr_0->m_seq_ids.Insert(7);
      EXPECT_EQ(1, ptr_0->m_seq_ids.GetSize());
      EXPECT_TRUE(ptr_0->m_seq_ids.Has(7));
      id = ptr_0->ID();
    }
    {
      cLabel* ptr_1 = idx.Get(id);
      EXPECT_EQ(1, ptr_1->m_seq_ids.GetSize());
      EXPECT_TRUE(ptr_1->m_seq_ids.Has(7));
    }
  }

  TEST(cLabelIdx, delete_by_id){
    cLabelIdx idx;
    /* Initial index size should be zero. */
    EXPECT_EQ(0, idx.GetSize());
    /* Delete by nonsense ID should fail. */
    EXPECT_FALSE(idx.Delete(-1));
    /* Delete by bad ID should fail. */
    EXPECT_FALSE(idx.Delete(0));
    /* Index size should still be zero. */
    EXPECT_EQ(0, idx.GetSize());
    /* Five should be ID of first object. */
    EXPECT_EQ(5, idx.Insert(5)->ID());
    /* Index size should now be one. */
    EXPECT_EQ(1, idx.GetSize());
    /* First delete by ID 5 should succeed. */
    EXPECT_TRUE(idx.Delete(5));
    /* After delete of sole object, index size should be zero again. */
    EXPECT_EQ(0, idx.GetSize());
    /* Second delete by ID 0 should fail, since object was already deleted. */
    EXPECT_FALSE(idx.Delete(5));
    /* Index size should still be zero. */
    EXPECT_EQ(0, idx.GetSize());
  }

  TEST(cLabelIdx, iterator){
    cLabelIdx idx;
    int i = 0;
    /* Create and verify three objects in index. */
    EXPECT_EQ(2, idx.Insert(2)->ID());
    EXPECT_EQ(3, idx.Insert(3)->ID());
    EXPECT_EQ(5, idx.Insert(5)->ID());
    EXPECT_EQ(3, idx.GetSize());
    /* Should iterate over three IDs. */
    i = 0;
    for (cLabelIdx::Iterator it = idx.Begin(); it.Next(); ){
      i+=1;
      EXPECT_TRUE(it.Get());
      EXPECT_LE(2, it.ID());
      EXPECT_GE(5, it.ID());
    }
    EXPECT_EQ(3, i);
    /* Delete and verify one object from index. */
    EXPECT_TRUE(idx.Delete(3));
    EXPECT_EQ(2, idx.GetSize());
    EXPECT_FALSE(idx.Delete(3));
    EXPECT_EQ(2, idx.GetSize());
    /* Should iterate over two IDs. */
    i = 0;
    for (cLabelIdx::Iterator it = idx.Begin(); it.Next(); ){
      i+=1;
      EXPECT_TRUE(it.Get());
      EXPECT_LE(2, it.ID());
      EXPECT_NE(3, it.ID());
      EXPECT_GE(5, it.ID());
    }
    EXPECT_EQ(2, i);
    /* Regression: Get with bad ID should not accidentally add to index! */
    /* Should iterate over two IDs. */
    EXPECT_FALSE(idx.Get(-1));
    EXPECT_EQ(2, idx.GetSize());
    EXPECT_FALSE(idx.Get(1));
    EXPECT_EQ(2, idx.GetSize());
    EXPECT_FALSE(idx.Get(3));
    EXPECT_EQ(2, idx.GetSize());
    i = 0;
    for (cLabelIdx::Iterator it = idx.Begin(); it.Next(); ){
      i+=1;
      EXPECT_TRUE(it.Get());
      EXPECT_LE(2, it.ID());
      EXPECT_NE(3, it.ID());
      EXPECT_GE(5, it.ID());
    }
    EXPECT_EQ(2, i);
  }
}


namespace nFSMDBTests {

  class cFSMDBTestFixture : public cFSMDB {
  public:
    int InsertSequence(const Apto::String &sequence) { return cFSMDB::InsertSequence(sequence); }
    void UnlinkSeqLbls(int seq_id) { cFSMDB::UnlinkSeqLbls(seq_id); }
    bool RemoveSequence(int seq_id) { return cFSMDB::RemoveSequence(seq_id); }
    int CreateStrand(const Apto::String &sequence) { return cFSMDB::CreateStrand(sequence); }
  };


//  TEST(cFSMdb, instantiation){
//    cFSMDBTestFixture db;
//    /* Simple check to make sure instantiation works. */
//    /* Call methods to ensure instantiation of template class. */
//    /* Indices should initially be empty, so have zero size! */
//    EXPECT_EQ(0, db.m_seq_idx.GetSize());
//    EXPECT_EQ(0, db.m_strand_idx.GetSize());
//    EXPECT_EQ(0, db.m_label_idx.GetSize());
//  }


//  TEST(cFSMdb, brainstorm){
//    cFSMDBTestFixture db;
//    int seq_id = db.m_seq_idx.Insert("abcdxyzabcdxyzabcd");
//    /*
//    Scenario:
//    - strand-to-seq relation is many-to-one.
//    - seq-to-label relation is many-to-many.
//    - If we delete a strand, and it is the last strand associated with a
//      sequence, we should also delete the sequence.
//    - If we delete a sequence, and certain labels are associated with only
//      that sequence, we should also delete these labels.
//
//    - A cell/hardware/organism owns strands.
//      - FSMs ask the cell/hardware/organism to delete or create strands.
//      - FSMs can modify strands.
//    - A cell/hardware/organism owns FSMs.
//      - FSMs ask the cell/hardware/organism to delete or create FSMs.
//      - FSMs can modify FSMs.
//    - A strand owns one sequence; but many strands may share a sequence.
//    - An FSM owns one coding; but many FMSs may share a coding.
//    - Sequences and FSM codings own labels; many sequences and FSM codings
//      may share many labels.
//    - Strands and FSMs own bindings. A strand or FSM can own many bindings. A
//      binding is owned by two objects, which may be strands or FSMs.
//    */
//    //cStrand.
//  }


//  TEST(cFSMdb, sequence_accessors){
//    cFSMDBTestFixture db;
//    Apto::String seq("abcdxyzabcdxyzabcd");
//  
//    //EXPECT_FALSE(db.HasSequence(-1));
//    //EXPECT_FALSE(db.HasSequence(0));
//    //EXPECT_FALSE(db.HasSequence(Apto::SmartPtr<cSequence>());
//    //EXPECT_FALSE(db.HasSequence(seq));
//  
//    //EXPECT_EQ("", db.GetSequenceString(-1));
//    //EXPECT_EQ("", db.GetSequenceString(0));
//    //EXPECT_EQ("", db.GetSequenceString(Apto::SmartPtr<cSequence>()));
//    //EXPECT_EQ((cLabel*)0, db.GetSequence(-1));
//    //EXPECT_EQ((cLabel*)0, db.GetSequence(0));
//    //EXPECT_EQ(-1, db.GetSequenceID(Apto::SmartPtr<cSequence>());
//    //EXPECT_EQ(-1, db.GetSequenceID(seq);
//  
//    //EXPECT_EQ(-1, db.SequenceRefCt(-1));
//    //EXPECT_EQ(-1, db.SequenceRefCt(0));
//    //EXPECT_EQ(-1, db.SequenceRefCt(Apto::SmartPtr<cSequence>());
//    //EXPECT_EQ(-1, db.SequenceRefCt(seq));
//  
//    //int id = db.InsertSequence(seq);
//    //Apto::SmartPtr<cSequence> ptr = db.GetSequence(id);
//  
//    //EXPECT_TRUE(db.HasSequence(id));
//    //EXPECT_TRUE(db.HasSequence(ptr));
//    //EXPECT_TRUE(db.HasSequence(seq));
//  
//    //EXPECT_EQ(seq, db.GetSequenceString(id));
//    //EXPECT_EQ(seq, db.GetSequenceString(ptr));
//    //EXPECT_EQ(ptr, db.GetSequence(id));
//    //EXPECT_EQ(ptr, db.GetSequence(seq));
//    //EXPECT_EQ(id, db.GetSequenceID(ptr));
//    //EXPECT_EQ(id, db.GetSequenceID(seq));
//  
//    //EXPECT_EQ(0, db.SequenceRefCt(id));
//    //EXPECT_EQ(0, db.SequenceRefCt(ptr));
//    //EXPECT_EQ(0, db.SequenceRefCt(seq));
//  }


//  TEST(cFSMdb, label_accessors){
//    cFSMDBTestFixture db;
//  
//    EXPECT_FALSE(db.HasLabel(5));
//    EXPECT_FALSE(db.HasLabel((cLabel*)0));
//  
//    EXPECT_EQ((cLabel*)0, db.GetLabel(5));
//    EXPECT_EQ(-1, db.GetLabelID((cLabel*)0));
//  
//    EXPECT_EQ(-1, db.LabelRefCt(5));
//    EXPECT_EQ(-1, db.LabelRefCt((cLabel*)0));
//  
//    int lbl_id = 300;
//    int lbl_pos = 7;
//    cLabel* lbl_ptr = db.InsertLabel(lbl_id);
//    /* Sanity check. */
//    EXPECT_NE((cLabel*)0, lbl_ptr);
//  
//    EXPECT_TRUE(db.HasLabel(lbl_id));
//    EXPECT_TRUE(db.HasLabel(lbl_ptr));
//  
//    EXPECT_EQ(lbl_ptr, db.GetLabel(lbl_id));
//    EXPECT_EQ(lbl_id, db.GetLabelID(lbl_ptr));
//  
//    EXPECT_EQ(0, db.LabelRefCt(lbl_id));
//    EXPECT_EQ(0, db.LabelRefCt(lbl_ptr));
//  
//    Apto::String seq_0("aaabaacaadabbabcabdacbaccacdadbadcaddbbbcbbdbccbcdbdcbddcccdcddd");
//    Apto::String seq_1("abcdxyzabcdxyzabcd");
//    Apto::String seq_2("abcdabcdabcd");
//    Apto::String seq_3("abcdxyzabcd");
//    Apto::String seq_4("abcdabcd");
//    int seq_id_0 = db.InsertSequence(seq_0);
//    int seq_id_1 = db.InsertSequence(seq_1);
//    int seq_id_2 = db.InsertSequence(seq_2);
//    int seq_id_3 = db.InsertSequence(seq_3);
//    int seq_id_4 = db.InsertSequence(seq_4);
//    db.LinkSeqLblPos(seq_id_0, lbl_id, lbl_pos);
//    db.LinkSeqLblPos(seq_id_1, lbl_id, lbl_pos);
//    db.LinkSeqLblPos(seq_id_2, lbl_id, lbl_pos);
//    db.LinkSeqLblPos(seq_id_3, lbl_id, lbl_pos);
//    db.LinkSeqLblPos(seq_id_4, lbl_id, lbl_pos);
//  
//    EXPECT_EQ(5, db.LabelRefCt(lbl_id));
//    EXPECT_EQ(5, db.LabelRefCt(lbl_ptr));
//  
//    /* Unlink nonsense position; refct shouldn't change. */
//    db.UnlinkSeqLblPos(seq_id_4, lbl_id, -1);
//    EXPECT_EQ(5, db.LabelRefCt(lbl_id));
//    EXPECT_EQ(5, db.LabelRefCt(lbl_ptr));
//  
//    db.UnlinkSeqLblPos(seq_id_4, lbl_id, lbl_pos);
//    EXPECT_EQ(4, db.LabelRefCt(lbl_id));
//    EXPECT_EQ(4, db.LabelRefCt(lbl_ptr));
//  
//    db.UnlinkSeqLbl(seq_id_3, lbl_id);
//    EXPECT_EQ(3, db.LabelRefCt(lbl_id));
//    EXPECT_EQ(3, db.LabelRefCt(lbl_ptr));
//  
//    db.UnlinkSeqLbls(seq_id_2);
//    EXPECT_EQ(2, db.LabelRefCt(lbl_id));
//    EXPECT_EQ(2, db.LabelRefCt(lbl_ptr));
//  
//    db.UnlinkSeqLbls(seq_id_1);
//    EXPECT_EQ(1, db.LabelRefCt(lbl_id));
//    EXPECT_EQ(1, db.LabelRefCt(lbl_ptr));
//  
//    db.RemoveSequence(seq_id_0);
//    EXPECT_FALSE(db.HasLabel(lbl_id));
//    EXPECT_FALSE(db.HasLabel(lbl_ptr));
//    EXPECT_EQ(-1, db.LabelRefCt(lbl_id));
//    EXPECT_EQ(-1, db.LabelRefCt(lbl_ptr));
//  
//    EXPECT_FALSE(db.HasSequence(seq_id_0));
//  }


//  TEST(cLabel, label_ref_ct) {
//    cFSMDBTestFixture db;
//    cLabel* lbl_ptr = db.InsertLabel(300);
//    EXPECT_EQ(0, lbl_ptr->GetSeqCount());
//    lbl_ptr->InsertSeq(0);
//    EXPECT_EQ(1, lbl_ptr->GetSeqCount());
//    lbl_ptr->InsertSeq(0);
//    EXPECT_EQ(1, lbl_ptr->GetSeqCount());
//    lbl_ptr->InsertSeq(1);
//    EXPECT_EQ(2, lbl_ptr->GetSeqCount());
//    lbl_ptr->RemoveSeq(2);
//    EXPECT_EQ(2, lbl_ptr->GetSeqCount());
//    lbl_ptr->RemoveSeq(1);
//    EXPECT_EQ(1, lbl_ptr->GetSeqCount());
//    lbl_ptr->RemoveSeq(1);
//    EXPECT_EQ(1, lbl_ptr->GetSeqCount());
//    lbl_ptr->RemoveSeq(0);
//    EXPECT_EQ(0, lbl_ptr->GetSeqCount());
//  }


//  TEST(cFSMdb, sequence_is_referenced){
//  }



  TEST(cFSMDB, brainstorm_0){
    cFSMDBTestFixture db;
    /* Create first seq. */
    Apto::String seq_0("abcdxyzabcdxyzabcd");
    ///* Sanity: this should be an unseen sequence. */
    EXPECT_FALSE(db.m_seqs.Has(seq_0));
    /*
    Add the sequence; since it hasn't been seen, this should scan and index the
    sequence for labels.
    */
    int seq_id_0 = db.InsertSequence(seq_0);
    EXPECT_EQ("abcdxyzabcdxyzabcd", db.m_seqs.GetString(seq_id_0));
    /* Sanity: it should now be a seen sequence. */
    EXPECT_TRUE(db.m_seqs.Has(seq_0));
    cSequence* seq_ptr_0 = db.m_seqs.Get(seq_id_0);
    /* Sanity: make sure seq_ptr_0 is valid. */
    EXPECT_TRUE(seq_ptr_0);
    /*
    Adding the sequence a second time should not result in reindexing, which is
    reflected in the verification below of label info.
    */
    EXPECT_EQ(seq_id_0, db.InsertSequence(seq_0));
  
    /* Verify stored label info in sequence and label objects. */
    EXPECT_EQ(9, seq_ptr_0->m_lbl_sites.GetSize());
    int lbl_id = -1;
    int lbl_pos = -1;
    //Apto::SmartPtr<cLabel> lbl_ptr;
    cLabel* lbl_ptr;
    //Apto::Map<int, Apto::Set<int> >::KeyIterator lbl_it = seq_ptr_0->m_lbl_sites.Keys();  
    Apto::Map<int, Apto::Array<int> >::KeyIterator lbl_it = seq_ptr_0->m_lbl_sites.Keys();  
  
    lbl_it.Next(); lbl_id = *lbl_it.Get(); lbl_ptr = db.m_lbls.Get(lbl_id); 
    EXPECT_EQ(db.m_label_utils.Seq2ID("a"), lbl_id);
    EXPECT_EQ(3, seq_ptr_0->m_lbl_sites[lbl_id].GetSize());
    EXPECT_EQ(1, lbl_ptr->m_seq_ids.GetSize());
    EXPECT_TRUE(lbl_ptr->m_seq_ids.Has(seq_id_0));
    //Apto::Set<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin(); {
    {
      Apto::Array<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin();
      pos_it.Next(); EXPECT_EQ(0, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(7, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(14, *pos_it.Get());
    }
    lbl_it.Next(); lbl_id = *lbl_it.Get(); lbl_ptr = db.m_lbls.Get(lbl_id); 
    EXPECT_EQ(db.m_label_utils.Seq2ID("b"), lbl_id);
    EXPECT_EQ(3, seq_ptr_0->m_lbl_sites[lbl_id].GetSize());
    EXPECT_EQ(1, lbl_ptr->m_seq_ids.GetSize());
    EXPECT_TRUE(lbl_ptr->m_seq_ids.Has(seq_id_0));
    {
      Apto::Array<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin();
      pos_it.Next(); EXPECT_EQ(1, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(8, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(15, *pos_it.Get());
    }
    lbl_it.Next(); lbl_id = *lbl_it.Get(); lbl_ptr = db.m_lbls.Get(lbl_id); 
    EXPECT_EQ(db.m_label_utils.Seq2ID("bcd"), lbl_id);
    EXPECT_EQ(3, seq_ptr_0->m_lbl_sites[lbl_id].GetSize());
    EXPECT_EQ(1, lbl_ptr->m_seq_ids.GetSize());
    EXPECT_TRUE(lbl_ptr->m_seq_ids.Has(seq_id_0));
    {
      Apto::Array<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin();
      pos_it.Next(); EXPECT_EQ(1, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(8, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(15, *pos_it.Get());
    }
    lbl_it.Next(); lbl_id = *lbl_it.Get(); lbl_ptr = db.m_lbls.Get(lbl_id); 
    EXPECT_EQ(db.m_label_utils.Seq2ID("c"), lbl_id);
    EXPECT_EQ(3, seq_ptr_0->m_lbl_sites[lbl_id].GetSize());
    EXPECT_EQ(1, lbl_ptr->m_seq_ids.GetSize());
    EXPECT_TRUE(lbl_ptr->m_seq_ids.Has(seq_id_0));
    {
      Apto::Array<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin();
      pos_it.Next(); EXPECT_EQ(2, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(9, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(16, *pos_it.Get());
    }
    lbl_it.Next(); lbl_id = *lbl_it.Get(); lbl_ptr = db.m_lbls.Get(lbl_id); 
    EXPECT_EQ(db.m_label_utils.Seq2ID("abc"), lbl_id);
    EXPECT_EQ(3, seq_ptr_0->m_lbl_sites[lbl_id].GetSize());
    EXPECT_EQ(1, lbl_ptr->m_seq_ids.GetSize());
    EXPECT_TRUE(lbl_ptr->m_seq_ids.Has(seq_id_0));
    {
      Apto::Array<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin();
      pos_it.Next(); EXPECT_EQ(0, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(7, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(14, *pos_it.Get());
    }
    lbl_it.Next(); lbl_id = *lbl_it.Get(); lbl_ptr = db.m_lbls.Get(lbl_id); 
    EXPECT_EQ(db.m_label_utils.Seq2ID("d"), lbl_id);
    EXPECT_EQ(3, seq_ptr_0->m_lbl_sites[lbl_id].GetSize());
    EXPECT_EQ(1, lbl_ptr->m_seq_ids.GetSize());
    EXPECT_TRUE(lbl_ptr->m_seq_ids.Has(seq_id_0));
    {
      Apto::Array<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin();
      pos_it.Next(); EXPECT_EQ(3, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(10, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(17, *pos_it.Get());
    }
    lbl_it.Next(); lbl_id = *lbl_it.Get(); lbl_ptr = db.m_lbls.Get(lbl_id); 
    EXPECT_EQ(db.m_label_utils.Seq2ID("ab"), lbl_id);
    EXPECT_EQ(3, seq_ptr_0->m_lbl_sites[lbl_id].GetSize());
    EXPECT_EQ(1, lbl_ptr->m_seq_ids.GetSize());
    EXPECT_TRUE(lbl_ptr->m_seq_ids.Has(seq_id_0));
    {
      Apto::Array<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin();
      pos_it.Next(); EXPECT_EQ(0, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(7, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(14, *pos_it.Get());
    }
    lbl_it.Next(); lbl_id = *lbl_it.Get(); lbl_ptr = db.m_lbls.Get(lbl_id); 
    EXPECT_EQ(db.m_label_utils.Seq2ID("bc"), lbl_id);
    EXPECT_EQ(3, seq_ptr_0->m_lbl_sites[lbl_id].GetSize());
    EXPECT_EQ(1, lbl_ptr->m_seq_ids.GetSize());
    EXPECT_TRUE(lbl_ptr->m_seq_ids.Has(seq_id_0));
    {
      Apto::Array<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin();
      pos_it.Next(); EXPECT_EQ(1, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(8, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(15, *pos_it.Get());
    }
    lbl_it.Next(); lbl_id = *lbl_it.Get(); lbl_ptr = db.m_lbls.Get(lbl_id); 
    EXPECT_EQ(db.m_label_utils.Seq2ID("cd"), lbl_id);
    EXPECT_EQ(3, seq_ptr_0->m_lbl_sites[lbl_id].GetSize());
    EXPECT_EQ(1, lbl_ptr->m_seq_ids.GetSize());
    EXPECT_TRUE(lbl_ptr->m_seq_ids.Has(seq_id_0));
    {
      Apto::Array<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin();
      pos_it.Next(); EXPECT_EQ(2, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(9, *pos_it.Get());
      pos_it.Next(); EXPECT_EQ(16, *pos_it.Get());
    }
  
    /* To help regenerate the above: */
    //for (Apto::Map<int, Apto::Array<int> >::KeyIterator lbl_it = seq_ptr_0->m_lbl_sites.Keys(); lbl_it.Next();) {
    //  int lbl_id = *lbl_it.Get();
    //  cout << "lbl_it.Next(); lbl_id = *lbl_it.Get(); lbl_ptr = db.m_lbls.Get(lbl_id);" << endl;
    //  cout << "EXPECT_EQ(db.m_label_utils.Seq2ID(\"" << db.m_label_utils.ID2Seq(lbl_id) << "\"), lbl_id);" << endl;
    //  cout << "EXPECT_EQ(" << seq_ptr_0->m_lbl_sites[lbl_id].GetSize() << ", seq_ptr_0->m_lbl_sites[lbl_id].GetSize());" << endl;
    //  cout << "EXPECT_EQ(" << lbl_ptr->m_seq_ids.GetSize() << ", lbl_ptr->m_seq_ids.GetSize());" << endl;
    //  cout << "EXPECT_TRUE(lbl_ptr->m_seq_ids.Has(seq_id_0));" << endl;
    //  cout << "{" << endl;
    //  cout << "  Apto::Array<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin(); {" << endl;
    //  for (Apto::Array<int>::Iterator pos_it = seq_ptr_0->m_lbl_sites[lbl_id].Begin(); pos_it.Next();) {
    //    int lbl_pos = *pos_it.Get();
    //    cout << "  pos_it.Next(); EXPECT_EQ(" << lbl_pos << ", *pos_it.Get());" << endl;
    //  }
    //  cout << "}" << endl;
    //}
  
    /*
    Need to lookup:
    - Given label ID, list of matching sequences and positions
      // Array<bSeqHit> *hits = db.GetSequenceHits(lbl_id)
      - And subsequently from sequences to strands
      - Why? Because the process of binding a strand locates strands via labels.
    - Given label ID, list of matching FSM models and positions
      // Array<cFSMHit> *hits = db.getFSMHits(lbl_id)
      - And subsequently from FSM models to FSMs
      - Why? Because the process of binding an FSM locates n FSMs via labels.
    - Given a sequence, list of matching labels.
      - Why? Because when the last sequence or FSM model referencing a label is
        deleted, the label needs to be deleted too.
        - Array<bLblHit> *hits = db.getLblHits(seq_id)
        - for (int i=0; i < hits->GetSize(); i++) {
            
          }
    - Given a sequence, list of matching strands.
      - Why? To locate a strand for binding.
    - Given a strand, corresponding sequence.
      - Why? When last strand reference a sequence is deleted, the sequence needs
        to be deleted too.
    - Given an FSM model, list of matching FSMs.
      - Why? To locate an FSM for binding.
    - Given an FSM, corresponding FSM model.
      - Why? When last FSM referencing an FSM model is deleted, the FSM model
        needs to be deleted too.
    */
    /*
    - Option:
      - Map<lbl_id, Array<seq_ids> >
      - Map<seq_ids, Array<Pair<lbl_id, lbl_pos> >
      - This won't work, because I need to be able to remove the entry for a
        seq_ids from the first map.
    - Option:
      - Map<lbl_id, Map<seq_ids, Array<lbl_pos> > >
      - Map<seq_ids, Array<Pair<lbl_id, lbl_pos> >
      - This would work, but it has redundant info. How much can I limit redundancy?
    - Option:
      - Map<lbl_id, Map<seq_ids, Array<lbl_pos> > >
      - Map<seq_ids, Map<lbl_id, Array<lbl_pos> > >
      - This would work, but it has redundant info. How much can I limit redundancy?
    - Option:
      - Map<lbl_id, Set<seq_ids> >
      - Within sequence: Map<lbl_ids, Array<lbl_pos> >
      - This would work.
    - Option:
      - Within lable: Set<seq_ids>
      - Within sequence: Map<lbl_ids, Array<lbl_pos> >
      - This would work, and would use less memory.
    */
    ///* Create first strand. */
    //cStrand* strand_ptr_0 = db.m_strand_idx.Create();
    //int strand_id_0 = db.m_strand_idx.GetID(strand_ptr_0);
    ///* Associate with exactly one seq. */
    //strand_ptr_0->m_seq_id = seq_id_0;
    ///* Create second strand. */
    //cStrand* strand_ptr_1 = db.m_strand_idx.Create();
    //int strand_id_1 = db.m_strand_idx.GetID(strand_ptr_1);
    ///* Associate with exactly one seq. */
    //strand_ptr_1->m_seq_id = seq_id_0;
  
    /* Given sequence, lookup strands. */
    /* Given label, lookup strands. */
  }


  TEST(cFSMDB, brainstorm_1) {
    cFSMDBTestFixture db;
    /*
    This is a debruijn sequence containing:
    - one each of all triplet labels with 'a', 'b', 'c', or 'd';
    - four each of all couplet labels with 'a', 'b', 'c', or 'd';
    - 64 each of all singlet labels with 'a', 'b', 'c', or 'd';
    - one extra 'aa' and two extra 'a's (because I made the tail of the sequence
      overlap with the head by two characters, in order to ensure that every
      triplet is present).
    */
    Apto::String seq_0("aaabaacaadabbabcabdacbaccacdadbadcaddbbbcbbdbccbcdbdcbddcccdcdddaa");
    int seq_id_0 = db.InsertSequence(seq_0);
    EXPECT_EQ(1, db.m_seqs.GetSize());
    /*
    There are 64 unique labels of length three; 16 of length two; and four of
    length one; so there are a total of 84 unique labels appearing in the
    debruijn sequence.
    */
    EXPECT_EQ(84, db.m_lbls.GetSize());
    db.UnlinkSeqLbls(seq_id_0);
    EXPECT_EQ(0, db.m_lbls.GetSize());
  }


  TEST(cFSMDB, brainstorm_2) {
    cFSMDBTestFixture db;
    /*
    This is a debruijn sequence containing:
    - one each of all triplet labels with 'a', 'b', 'c', or 'd';
    - four each of all couplet labels with 'a', 'b', 'c', or 'd';
    - 64 each of all singlet labels with 'a', 'b', 'c', or 'd';
    - one extra 'aa' and two extra 'a's (because I made the tail of the sequence
      overlap with the head by two characters, in order to ensure that every
      triplet is present).
    */
    Apto::String seq_0("aaabaacaadabbabcabdacbaccacdadbadcaddbbbcbbdbccbcdbdcbddcccdcdddaa");
    int strand_id_0 = db.CreateStrand(seq_0);
    EXPECT_EQ(1, db.m_bindables.GetSize());
    EXPECT_EQ(1, db.m_seqs.GetSize());
    /*
    There are 64 unique labels of length three; 16 of length two; and four of
    length one; so there are a total of 84 unique labels appearing in the
    debruijn sequence.
    */
    EXPECT_EQ(84, db.m_lbls.GetSize());
    int strand_id_1 = db.CreateStrand(seq_0);
    EXPECT_EQ(2, db.m_bindables.GetSize());
    EXPECT_EQ(1, db.m_seqs.GetSize());
    EXPECT_EQ(84, db.m_lbls.GetSize());
    db.RemoveStrand(strand_id_1);
    EXPECT_EQ(1, db.m_bindables.GetSize());
    EXPECT_EQ(1, db.m_seqs.GetSize());
    EXPECT_EQ(84, db.m_lbls.GetSize());
    db.RemoveStrand(strand_id_0);
    EXPECT_EQ(0, db.m_bindables.GetSize());
    EXPECT_EQ(0, db.m_seqs.GetSize());
    EXPECT_EQ(0, db.m_lbls.GetSize());
  }


  TEST(cFSMDB, brainstorm_3) {
    cFSMDBTestFixture db;

  }
}


namespace nBasicLabelUtilsTests {
  TEST(cBasicLabelUtils, Rotate) {
    cBasicLabelUtils u;
    int id = u.Seq2ID("abcd");
    int new_id = u.Rotate(id);
    EXPECT_EQ(u.Seq2ID("abcd"), id);
    EXPECT_EQ(u.Seq2ID("cdab"), new_id);
  }

  TEST(cBasicLabelUtils, ReverseComplement) {
    cBasicLabelUtils u;
    int id = u.Seq2ID("abcd");
    int new_id = u.ReverseComplement(id);
    EXPECT_EQ(u.Seq2ID("abcd"), id);
    EXPECT_EQ(u.Seq2ID("badc"), new_id);
  }
}


namespace nAptoSchedulerDynamicTests {
  class ProbabilisticDynamicTestFixture : public Apto::Scheduler::ProbabilisticDynamic {
  public:
    ProbabilisticDynamicTestFixture(int num_entries, Apto::SmartPtr<Apto::RNG::AvidaRNG> rng)
    : Apto::Scheduler::ProbabilisticDynamic(num_entries, rng)
    {}
    double ItemWeight(int entry_id) { return m_index.m_item_weight[entry_id]; }
    double SubtreeWeight(int entry_id) { return m_index.m_subtree_weight[entry_id]; }
  };

  TEST(ProbabilisticDynamic, brainstorm_0) {
    Apto::SmartPtr<Apto::RNG::AvidaRNG> rng(new Apto::RNG::AvidaRNG);
    ProbabilisticDynamicTestFixture scheduler(5, rng);

    for (int i = 0; i < scheduler.GetSize(); i++) { scheduler.AdjustPriority(i, i); }
    EXPECT_EQ(0, scheduler.ItemWeight(0)); EXPECT_EQ(10, scheduler.SubtreeWeight(0));
    EXPECT_EQ(1, scheduler.ItemWeight(1)); EXPECT_EQ(8, scheduler.SubtreeWeight(1));
    EXPECT_EQ(2, scheduler.ItemWeight(2)); EXPECT_EQ(2, scheduler.SubtreeWeight(2));
    EXPECT_EQ(3, scheduler.ItemWeight(3)); EXPECT_EQ(3, scheduler.SubtreeWeight(3));
    EXPECT_EQ(4, scheduler.ItemWeight(4)); EXPECT_EQ(4, scheduler.SubtreeWeight(4));

    /* To generate above EXPECTs: */
    //for (int i = 0; i < scheduler.GetSize(); i++) {
    //  cout << "EXPECT_EQ(" << scheduler.ItemWeight(i) << ", scheduler.ItemWeight(" << i << ")); EXPECT_EQ(" << scheduler.SubtreeWeight(i) << ", scheduler.SubtreeWeight(" << i << "));" << endl;
    //}
  }

  TEST(ProbabilisticDynamic, brainstorm_1) {
    Apto::SmartPtr<Apto::RNG::AvidaRNG> rng(new Apto::RNG::AvidaRNG);
    ProbabilisticDynamicTestFixture scheduler(20, rng);

    for (int i = 0; i < scheduler.GetSize(); i++) { scheduler.AdjustPriority(i, i); }
    EXPECT_EQ(0, scheduler.ItemWeight(0)); EXPECT_EQ(190, scheduler.SubtreeWeight(0));
    EXPECT_EQ(1, scheduler.ItemWeight(1)); EXPECT_EQ(127, scheduler.SubtreeWeight(1));
    EXPECT_EQ(2, scheduler.ItemWeight(2)); EXPECT_EQ(63, scheduler.SubtreeWeight(2));
    EXPECT_EQ(3, scheduler.ItemWeight(3)); EXPECT_EQ(84, scheduler.SubtreeWeight(3));
    EXPECT_EQ(4, scheduler.ItemWeight(4)); EXPECT_EQ(42, scheduler.SubtreeWeight(4));
    EXPECT_EQ(5, scheduler.ItemWeight(5)); EXPECT_EQ(28, scheduler.SubtreeWeight(5));
    EXPECT_EQ(6, scheduler.ItemWeight(6)); EXPECT_EQ(33, scheduler.SubtreeWeight(6));
    EXPECT_EQ(7, scheduler.ItemWeight(7)); EXPECT_EQ(38, scheduler.SubtreeWeight(7));
    EXPECT_EQ(8, scheduler.ItemWeight(8)); EXPECT_EQ(43, scheduler.SubtreeWeight(8));
    EXPECT_EQ(9, scheduler.ItemWeight(9)); EXPECT_EQ(28, scheduler.SubtreeWeight(9));
    EXPECT_EQ(10, scheduler.ItemWeight(10)); EXPECT_EQ(10, scheduler.SubtreeWeight(10));
    EXPECT_EQ(11, scheduler.ItemWeight(11)); EXPECT_EQ(11, scheduler.SubtreeWeight(11));
    EXPECT_EQ(12, scheduler.ItemWeight(12)); EXPECT_EQ(12, scheduler.SubtreeWeight(12));
    EXPECT_EQ(13, scheduler.ItemWeight(13)); EXPECT_EQ(13, scheduler.SubtreeWeight(13));
    EXPECT_EQ(14, scheduler.ItemWeight(14)); EXPECT_EQ(14, scheduler.SubtreeWeight(14));
    EXPECT_EQ(15, scheduler.ItemWeight(15)); EXPECT_EQ(15, scheduler.SubtreeWeight(15));
    EXPECT_EQ(16, scheduler.ItemWeight(16)); EXPECT_EQ(16, scheduler.SubtreeWeight(16));
    EXPECT_EQ(17, scheduler.ItemWeight(17)); EXPECT_EQ(17, scheduler.SubtreeWeight(17));
    EXPECT_EQ(18, scheduler.ItemWeight(18)); EXPECT_EQ(18, scheduler.SubtreeWeight(18));
    EXPECT_EQ(19, scheduler.ItemWeight(19)); EXPECT_EQ(19, scheduler.SubtreeWeight(19));

    scheduler.Resize(5);
    EXPECT_EQ(0, scheduler.ItemWeight(0)); EXPECT_EQ(10, scheduler.SubtreeWeight(0));
    EXPECT_EQ(1, scheduler.ItemWeight(1)); EXPECT_EQ(8, scheduler.SubtreeWeight(1));
    EXPECT_EQ(2, scheduler.ItemWeight(2)); EXPECT_EQ(2, scheduler.SubtreeWeight(2));
    EXPECT_EQ(3, scheduler.ItemWeight(3)); EXPECT_EQ(3, scheduler.SubtreeWeight(3));
    EXPECT_EQ(4, scheduler.ItemWeight(4)); EXPECT_EQ(4, scheduler.SubtreeWeight(4));

    scheduler.Resize(20);
    for (int i = 5; i < scheduler.GetSize(); i++) { scheduler.AdjustPriority(i, i); }
    EXPECT_EQ(0, scheduler.ItemWeight(0)); EXPECT_EQ(190, scheduler.SubtreeWeight(0));
    EXPECT_EQ(1, scheduler.ItemWeight(1)); EXPECT_EQ(127, scheduler.SubtreeWeight(1));
    EXPECT_EQ(2, scheduler.ItemWeight(2)); EXPECT_EQ(63, scheduler.SubtreeWeight(2));
    EXPECT_EQ(3, scheduler.ItemWeight(3)); EXPECT_EQ(84, scheduler.SubtreeWeight(3));
    EXPECT_EQ(4, scheduler.ItemWeight(4)); EXPECT_EQ(42, scheduler.SubtreeWeight(4));
    EXPECT_EQ(5, scheduler.ItemWeight(5)); EXPECT_EQ(28, scheduler.SubtreeWeight(5));
    EXPECT_EQ(6, scheduler.ItemWeight(6)); EXPECT_EQ(33, scheduler.SubtreeWeight(6));
    EXPECT_EQ(7, scheduler.ItemWeight(7)); EXPECT_EQ(38, scheduler.SubtreeWeight(7));
    EXPECT_EQ(8, scheduler.ItemWeight(8)); EXPECT_EQ(43, scheduler.SubtreeWeight(8));
    EXPECT_EQ(9, scheduler.ItemWeight(9)); EXPECT_EQ(28, scheduler.SubtreeWeight(9));
    EXPECT_EQ(10, scheduler.ItemWeight(10)); EXPECT_EQ(10, scheduler.SubtreeWeight(10));
    EXPECT_EQ(11, scheduler.ItemWeight(11)); EXPECT_EQ(11, scheduler.SubtreeWeight(11));
    EXPECT_EQ(12, scheduler.ItemWeight(12)); EXPECT_EQ(12, scheduler.SubtreeWeight(12));
    EXPECT_EQ(13, scheduler.ItemWeight(13)); EXPECT_EQ(13, scheduler.SubtreeWeight(13));
    EXPECT_EQ(14, scheduler.ItemWeight(14)); EXPECT_EQ(14, scheduler.SubtreeWeight(14));
    EXPECT_EQ(15, scheduler.ItemWeight(15)); EXPECT_EQ(15, scheduler.SubtreeWeight(15));
    EXPECT_EQ(16, scheduler.ItemWeight(16)); EXPECT_EQ(16, scheduler.SubtreeWeight(16));
    EXPECT_EQ(17, scheduler.ItemWeight(17)); EXPECT_EQ(17, scheduler.SubtreeWeight(17));
    EXPECT_EQ(18, scheduler.ItemWeight(18)); EXPECT_EQ(18, scheduler.SubtreeWeight(18));
    EXPECT_EQ(19, scheduler.ItemWeight(19)); EXPECT_EQ(19, scheduler.SubtreeWeight(19));

    scheduler.Resize(5);
    EXPECT_EQ(0, scheduler.ItemWeight(0)); EXPECT_EQ(10, scheduler.SubtreeWeight(0));
    EXPECT_EQ(1, scheduler.ItemWeight(1)); EXPECT_EQ(8, scheduler.SubtreeWeight(1));
    EXPECT_EQ(2, scheduler.ItemWeight(2)); EXPECT_EQ(2, scheduler.SubtreeWeight(2));
    EXPECT_EQ(3, scheduler.ItemWeight(3)); EXPECT_EQ(3, scheduler.SubtreeWeight(3));
    EXPECT_EQ(4, scheduler.ItemWeight(4)); EXPECT_EQ(4, scheduler.SubtreeWeight(4));

    /* To generate above EXPECTs: */
    //for (int i = 0; i < scheduler.GetSize(); i++) {
    //  cout << "EXPECT_EQ(" << scheduler.ItemWeight(i) << ", scheduler.ItemWeight(" << i << ")); EXPECT_EQ(" << scheduler.SubtreeWeight(i) << ", scheduler.SubtreeWeight(" << i << "));" << endl;
    //}
  }

  TEST(ProbabilisticDynamic, brainstorm_2) {
    /* Seed random number generator. */
    std::srand(0);
    cFSMDB db;
    db.m_label_utils.m_max_label_size = 4;
    /* debruijn sequence with all triplets of a-d. */
    Apto::String seq_0("aaabaacaadabbabcabdacbaccacdadbadcaddbbbcbbdbccbcdbdcbddcccdcdddaa");
    Apto::String seq_1("abcdxyzabcdxyzabcd");
    Apto::SmartPtr<Apto::RNG::AvidaRNG> rng(new Apto::RNG::AvidaRNG);

    /* Load strand molecules into scheduler. */
    for (int i=0; i<5; i++) { db.CreateStrand(seq_0); }
    for (int i=0; i<5; i++) { db.CreateStrand(seq_1); }

    /* Molecule collisions. */
    for (int i=0; i<20; i++) { db.SingleCollision(); }
  }

  TEST(IntegratedDynamic, brainstorm_0) {
    Apto::Scheduler::IntegratedDynamic scheduler(5);
    Apto::Array<int> hits(0);
    int calls_to_next = 0;

    hits.Resize(scheduler.GetSize());
    hits.SetAll(0);
    for (int i = 0; i < scheduler.GetSize(); i++) {
      int priority = i + 1;
      scheduler.AdjustPriority(i, priority);
    }
    calls_to_next = 0;
    for (int i = 0; i < scheduler.GetSize(); i++) {
      int priority = i + 1;
      calls_to_next += priority;
    }
    for (int i = 0; i < calls_to_next; i++) {
      int next = scheduler.Next();
      hits[next] = hits[next] + 1;
    }
    for (int i = 0; i < scheduler.GetSize(); i++) { EXPECT_EQ(i+1, hits[i]); }

    scheduler.Resize(10);
    hits.Resize(scheduler.GetSize());
    hits.SetAll(0);
    for (int i = 5; i < scheduler.GetSize(); i++) {
      int priority = i + 1;
      scheduler.AdjustPriority(i, priority);
    }
    calls_to_next = 0;
    for (int i = 0; i < scheduler.GetSize(); i++) {
      int priority = i + 1;
      calls_to_next += priority;
    }
    for (int i = 0; i < calls_to_next; i++) { int next = scheduler.Next(); }
    for (int i = 0; i < calls_to_next; i++) {
      int next = scheduler.Next();
      hits[next] = hits[next] + 1;
    }
    for (int i = 0; i < scheduler.GetSize(); i++) { EXPECT_EQ(i+1, hits[i]); }

    scheduler.Resize(5);
    hits.Resize(scheduler.GetSize());
    hits.SetAll(0);
    calls_to_next = 0;
    for (int i = 0; i < scheduler.GetSize(); i++) {
      int priority = i + 1;
      calls_to_next += priority;
    }
    for (int i = 0; i < calls_to_next; i++) {
      int next = scheduler.Next();
      hits[next] = hits[next] + 1;
    }
    for (int i = 0; i < scheduler.GetSize(); i++) { EXPECT_EQ(i+1, hits[i]); }
  }
}


namespace nObjDeleteMeIdxTests {
  class Feline : public ObjBase {
  public:
    static int s_ct;
    int m_an_instance_variable;
  public:
    Feline(): m_an_instance_variable(0) { s_ct++; }
    virtual ~Feline() { s_ct--; }
  };
  
  class Lion : public Feline {
  public:
    static int s_ct;
  public:
    Lion(){ s_ct++; }
    virtual ~Lion() { s_ct--; }
  };
  
  class Tiger : public Feline {
  public:
    static int s_ct;
  public:
    Tiger(){ s_ct++; }
    virtual ~Tiger() { s_ct--; }
  };
  
  int Feline::s_ct(0);
  int Lion::s_ct(0);
  int Tiger::s_ct(0);

  void PrintCounts() {
    cout << "Feline::s_ct: " << Feline::s_ct << endl;
    cout << "Lion::s_ct: " << Lion::s_ct << endl;
    cout << "Tiger::s_ct: " << Tiger::s_ct << endl;
  }

  
  class aClass {
  public:
    int m_variable;
  };
  TEST(AptoMap, DISABLED_int_SmartPtr_regression) {
    Apto::Map<int, Apto::SmartPtr<aClass> > id2object;
    for (int id = 0; id < 100; id++) {
      id2object.Set(id, Apto::SmartPtr<aClass>(new aClass()));
      bool ok = true;
      for (Apto::Map<int, Apto::SmartPtr<aClass> >::KeyIterator it = id2object.Keys(); it.Next();) {
        int test_id = *it.Get();
        if (!id2object.Has(test_id)) {
          cout << "id2object.Has(" << test_id << ")): " << id2object.Has(test_id) << endl;
          ok = false;
        } 
        if (!id2object.Get(test_id)) {
          cout << "id2object.Get(" << test_id << ")): " << id2object.Get(test_id) << endl;
          ok = false;
        } 
      }
      if (!ok) { cout << "problem after creating id " << id << endl; }
      EXPECT_TRUE(ok);
    }
    Apto::SmartPtr<aClass> obj_23(id2object[23]), obj_46(id2object[46]), obj_69(id2object[69]), obj_92(id2object[92]);
    {
      int id = 0;
      id2object.Remove(id);
      bool ok = true;
      for (Apto::Map<int, Apto::SmartPtr<aClass> >::KeyIterator it = id2object.Keys(); it.Next();) {
        int test_id = *it.Get();
        if (!id2object.Has(test_id)) {
          cout << "id2object.Has(" << test_id << ")): " << id2object.Has(test_id) << endl;
          ok = false;
        } 
        if (!id2object.Get(test_id)) {
          cout << "id2object.Get(" << test_id << ")): " << id2object.Get(test_id) << endl;
          ok = false;
        } 
      }
      if (!ok) { cout << "problem after deleting id " << id << endl; }
      EXPECT_TRUE(ok);
    }
    EXPECT_TRUE(obj_23);
    EXPECT_TRUE(obj_46);
    EXPECT_TRUE(obj_69);
    EXPECT_TRUE(obj_92);
    EXPECT_EQ(obj_23, id2object[23]);
    EXPECT_EQ(obj_46, id2object[46]);
    EXPECT_EQ(obj_69, id2object[69]);
    EXPECT_EQ(obj_92, id2object[92]);
    for (int id = 1; id < 100; id++) {
      id2object.Remove(id);
      bool ok = true;
      for (Apto::Map<int, Apto::SmartPtr<aClass> >::KeyIterator it = id2object.Keys(); it.Next();) {
        int test_id = *it.Get();
        if (!id2object.Has(test_id)) {
          cout << "id2object.Has(" << test_id << ")): " << id2object.Has(test_id) << endl;
          ok = false;
        } 
        if (!id2object.Get(test_id)) {
          cout << "id2object.Get(" << test_id << ")): " << id2object.Get(test_id) << endl;
          ok = false;
        } 
      }
      if (!ok) { cout << "problem after deleting id " << id << endl; }
      EXPECT_TRUE(ok);
    }
  }

  TEST(ObjIdx, instantiation) {
    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);

    ObjIdx<Feline> idx;

    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);
  }

  TEST(ObjIdx, initial_state) {
    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);

    ObjIdx<Feline> idx;
    EXPECT_EQ(0, idx.GetSize());
    EXPECT_EQ(0, idx.GetMaxCt());

    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);
  }

  TEST(ObjIdx, create_and_delete) {
    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);

    ObjIdx<Feline> idx;
    Feline* ptr = idx.Create();
    EXPECT_NE((Feline*)(0), ptr);
    EXPECT_EQ(1, idx.GetSize());
    EXPECT_EQ(1, idx.GetMaxCt());
    EXPECT_EQ(0, ptr->ID());

    EXPECT_EQ(1, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);

    EXPECT_TRUE(idx.Delete(0));
    EXPECT_EQ(0, idx.GetSize());
    EXPECT_EQ(1, idx.GetMaxCt());
    EXPECT_FALSE(idx.Delete(0));

    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);
  }

  TEST(ObjIdx, polymorphic_create_and_delete) {
    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);

    ObjIdx<Feline> idx;
    Lion* ptr = idx.Create<Lion>();
    EXPECT_NE((Lion*)(0), ptr);
    EXPECT_EQ(1, idx.GetSize());
    EXPECT_EQ(1, idx.GetMaxCt());
    EXPECT_EQ(0, ptr->ID());

    EXPECT_EQ(1, Feline::s_ct);
    EXPECT_EQ(1, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);

    EXPECT_TRUE(idx.Delete(0));
    EXPECT_EQ(0, idx.GetSize());
    EXPECT_EQ(1, idx.GetMaxCt());
    EXPECT_FALSE(idx.Delete(0));

    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);
  }

  TEST(ObjIdx, destructor) {
    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);

    {
      ObjIdx<Feline> idx;
      Lion* ptr = idx.Create<Lion>();
      EXPECT_NE((Lion*)(0), ptr);
      EXPECT_EQ(1, idx.GetSize());
      EXPECT_EQ(1, idx.GetMaxCt());
      EXPECT_EQ(0, ptr->ID());

      EXPECT_EQ(1, Feline::s_ct);
      EXPECT_EQ(1, Lion::s_ct);
      EXPECT_EQ(0, Tiger::s_ct);
    }

    /* We didn't explicitly Delete(); did destructor do it for us? */
    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);
  }

  TEST(ObjIdx, has) {
    ObjIdx<Feline> idx;
    EXPECT_FALSE(idx.Has(0));
    Feline* ptr = idx.Create();
    int id = ptr->ID();
    EXPECT_EQ(0, id);
    EXPECT_TRUE(idx.Has(0));
    EXPECT_TRUE(idx.Delete(0));
    EXPECT_FALSE(idx.Has(0));
  }

  TEST(ObjIdx, has_with_polymorphic) {
    ObjIdx<Feline> idx;
    EXPECT_FALSE(idx.Has(0));
    Lion* ptr = idx.Create<Lion>();
    int id = ptr->ID();
    EXPECT_EQ(0, id);
    EXPECT_TRUE(idx.Has(0));
    EXPECT_TRUE(idx.Delete(0));
    EXPECT_FALSE(idx.Has(0));
  }

  TEST(ObjIdx, get) {
    ObjIdx<Feline> idx;
    EXPECT_EQ((Feline*)(0), idx.Get(0));
    Feline* ptr = idx.Create();
    EXPECT_NE((Feline*)(0), ptr);
    int id = ptr->ID();
    EXPECT_EQ(0, id);
    EXPECT_EQ(ptr, idx.Get(0));
    EXPECT_TRUE(idx.Delete(0));
  }

  TEST(ObjIdx, polymorphic_get) {
    ObjIdx<Feline> idx;
    EXPECT_EQ((Feline*)(0), idx.Get(0));
    EXPECT_EQ((Feline*)(0), idx.Get(1));

    Feline* feline_ptr = idx.Create();
    EXPECT_NE((Feline*)(0), feline_ptr);
    int feline_id = feline_ptr->ID();
    EXPECT_EQ(0, feline_id);
    EXPECT_EQ(feline_ptr, idx.Get(0));
    /* dynamic cast of Feline to Lion should fail. */
    EXPECT_EQ((Lion*)(0), idx.Get<Lion>(0));

    Lion* lion_ptr = idx.Create<Lion>();
    EXPECT_NE((Lion*)(0), lion_ptr);
    int lion_id = lion_ptr->ID();
    EXPECT_EQ(1, lion_id);
    EXPECT_EQ(lion_ptr, idx.Get<Lion>(1));
    /* dynamic cast of Lion to Feline should succeed. */
    EXPECT_EQ((Feline*)(lion_ptr), idx.Get<Feline>(1));
    EXPECT_EQ((Feline*)(lion_ptr), idx.Get(1));

    Tiger* tiger_ptr = idx.Create<Tiger>();
    EXPECT_NE((Tiger*)(0), tiger_ptr);
    int tiger_id = tiger_ptr->ID();
    EXPECT_EQ(2, tiger_id);
    EXPECT_EQ(tiger_ptr, idx.Get<Tiger>(2));
    /* dynamic cast of Tiger to Feline should succeed. */
    EXPECT_EQ((Feline*)(tiger_ptr), idx.Get<Feline>(2));
    EXPECT_EQ((Feline*)(tiger_ptr), idx.Get(2));
    /* dynamic cast of Tiger to Lion should fail. */
    EXPECT_EQ((Lion*)(0), idx.Get<Lion>(2));

    EXPECT_TRUE(idx.Delete(2));
    EXPECT_TRUE(idx.Delete(1));
    EXPECT_TRUE(idx.Delete(0));
  }

  TEST(ObjIdx, handling_get_nonsense){
    ObjIdx<Feline> idx;
    EXPECT_EQ(0, idx.GetSize());
    /* Try to get using a nonsense ID; should return null ptr. */
    EXPECT_FALSE(idx.Get(-1));
    /* Try to get using a bad ID; should return null ptr. */
    EXPECT_FALSE(idx.Get(0));
    /* Get with bad ptr should not accidentally add to index! */
    EXPECT_EQ(0, idx.GetSize());
  }

  TEST(ObjIdx, persistence){
    ObjIdx<Feline> idx;
    /* Changes to a managed object should persist. */
    int id = -1;
    {
      Feline* ptr_0 = idx.Create();
      /* Initial value of instance variable should be zero. */
      EXPECT_EQ(0, ptr_0->m_an_instance_variable);
      ptr_0->m_an_instance_variable = 5;
      id = ptr_0->ID();
    }
    {
      Feline* ptr_1 = idx.Get(id);
      /*
      After value of instance variable was changed to 5, and object reaccessed,
      value should still be 5.
      */
      EXPECT_EQ(5, ptr_1->m_an_instance_variable);
    }
    EXPECT_TRUE(idx.Delete(id));
  }

  TEST(ObjIdx, iterator) {
    ObjIdx<Feline> idx;
    int i=0;
    for (i=0; i<100; i++){ EXPECT_EQ(i, idx.Create()->ID()); }
    EXPECT_EQ(100, idx.GetSize());
    EXPECT_EQ(100, idx.GetMaxCt());
    i=0;
    for (ObjIdx<Feline>::Iterator it = idx.Begin(); it.Next();){
      Feline* ptr = it.Get();
      EXPECT_NE((Feline*)(0), ptr);
      EXPECT_EQ(i, it.ID());
      EXPECT_EQ(i, ptr->ID());
      i++;
    }

    /* Verify iterator skips deleted objects. */
    EXPECT_TRUE(idx.Delete(0));
    EXPECT_TRUE(idx.Delete(1));
    i=2; /* We expect to start at id 2. */
    for (ObjIdx<Feline>::Iterator it = idx.Begin(); it.Next();){
      Feline* ptr = it.Get();
      EXPECT_NE((Feline*)(0), ptr);
      EXPECT_EQ(i, it.ID());
      EXPECT_EQ(i, ptr->ID());
      i++;
    }

    for (i=2; i<100; i++){ EXPECT_TRUE(idx.Delete(i)); }
    EXPECT_EQ(0, idx.GetSize());
    EXPECT_EQ(100, idx.GetMaxCt());
  }

  TEST(ObjIdx, polymorphic_iterator) {
    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);

    ObjIdx<Feline> idx;
    EXPECT_EQ(0, idx.Create<Feline>()->ID());
    EXPECT_EQ(1, idx.Create<Lion>()->ID());
    EXPECT_EQ(2, idx.Create<Tiger>()->ID());
    EXPECT_EQ(3, idx.Create<Feline>()->ID());
    EXPECT_EQ(4, idx.Create<Lion>()->ID());
    EXPECT_EQ(5, idx.Create<Tiger>()->ID());
    EXPECT_EQ(6, idx.Create<Feline>()->ID());
    EXPECT_EQ(7, idx.Create<Lion>()->ID());
    EXPECT_EQ(8, idx.Create<Tiger>()->ID());
    EXPECT_EQ(9, idx.Create<Feline>()->ID());

    EXPECT_EQ(10, Feline::s_ct);
    EXPECT_EQ(3, Lion::s_ct);
    EXPECT_EQ(3, Tiger::s_ct);

    /* Iterates through all objects, testing for Feline, Lion, or Tiger. */
    ObjIdx<Feline>::Iterator feline_it = idx.Begin();
    EXPECT_EQ(0, feline_it.Next<Feline>()->ID());
    EXPECT_EQ(0, feline_it.Get<Feline>()->ID());
    EXPECT_EQ((Lion*)(0), feline_it.Get<Lion>());
    EXPECT_EQ((Tiger*)(0), feline_it.Get<Tiger>());
    EXPECT_EQ(1, feline_it.Next<Feline>()->ID());
    EXPECT_EQ(1, feline_it.Get<Feline>()->ID());
    EXPECT_EQ(1, feline_it.Get<Lion>()->ID());
    EXPECT_EQ((Tiger*)(0), feline_it.Get<Tiger>());
    EXPECT_EQ(2, feline_it.Next<Feline>()->ID());
    EXPECT_EQ(2, feline_it.Get<Feline>()->ID());
    EXPECT_EQ((Lion*)(0), feline_it.Get<Lion>());
    EXPECT_EQ(2, feline_it.Get<Tiger>()->ID());
    EXPECT_EQ(3, feline_it.Next<Feline>()->ID());
    EXPECT_EQ(3, feline_it.Get<Feline>()->ID());
    EXPECT_EQ((Lion*)(0), feline_it.Get<Lion>());
    EXPECT_EQ((Tiger*)(0), feline_it.Get<Tiger>());
    EXPECT_EQ(4, feline_it.Next<Feline>()->ID());
    EXPECT_EQ(4, feline_it.Get<Feline>()->ID());
    EXPECT_EQ(4, feline_it.Get<Lion>()->ID());
    EXPECT_EQ((Tiger*)(0), feline_it.Get<Tiger>());
    EXPECT_EQ(5, feline_it.Next<Feline>()->ID());
    EXPECT_EQ(5, feline_it.Get<Feline>()->ID());
    EXPECT_EQ((Lion*)(0), feline_it.Get<Lion>());
    EXPECT_EQ(5, feline_it.Get<Tiger>()->ID());
    EXPECT_EQ(6, feline_it.Next<Feline>()->ID());
    EXPECT_EQ(6, feline_it.Get<Feline>()->ID());
    EXPECT_EQ((Lion*)(0), feline_it.Get<Lion>());
    EXPECT_EQ((Tiger*)(0), feline_it.Get<Tiger>());
    EXPECT_EQ(7, feline_it.Next<Feline>()->ID());
    EXPECT_EQ(7, feline_it.Get<Feline>()->ID());
    EXPECT_EQ(7, feline_it.Get<Lion>()->ID());
    EXPECT_EQ((Tiger*)(0), feline_it.Get<Tiger>());
    EXPECT_EQ(8, feline_it.Next<Feline>()->ID());
    EXPECT_EQ(8, feline_it.Get<Feline>()->ID());
    EXPECT_EQ((Lion*)(0), feline_it.Get<Lion>());
    EXPECT_EQ(8, feline_it.Get<Tiger>()->ID());
    EXPECT_EQ(9, feline_it.Next<Feline>()->ID());
    EXPECT_EQ(9, feline_it.Get<Feline>()->ID());
    EXPECT_EQ((Lion*)(0), feline_it.Get<Lion>());
    EXPECT_EQ((Tiger*)(0), feline_it.Get<Tiger>());
    EXPECT_EQ((Feline*)(0), feline_it.Next<Feline>());
    EXPECT_EQ((Feline*)(0), feline_it.Get<Feline>());

    /* Iterates through all Lion objects. */
    ObjIdx<Feline>::Iterator lion_it = idx.Begin();
    EXPECT_EQ(1, lion_it.Next<Lion>()->ID());
    EXPECT_EQ(1, lion_it.Get<Lion>()->ID());
    EXPECT_EQ(4, lion_it.Next<Lion>()->ID());
    EXPECT_EQ(4, lion_it.Get<Lion>()->ID());
    EXPECT_EQ(7, lion_it.Next<Lion>()->ID());
    EXPECT_EQ(7, lion_it.Get<Lion>()->ID());
    EXPECT_EQ((Lion*)(0), lion_it.Next<Lion>());
    EXPECT_EQ((Lion*)(0), lion_it.Get<Lion>());

    /* Iterates through all Tiger objects. */
    ObjIdx<Feline>::Iterator tiger_it = idx.Begin();
    EXPECT_EQ(2, tiger_it.Next<Tiger>()->ID());
    EXPECT_EQ(2, tiger_it.Get<Tiger>()->ID());
    EXPECT_EQ(5, tiger_it.Next<Tiger>()->ID());
    EXPECT_EQ(5, tiger_it.Get<Tiger>()->ID());
    EXPECT_EQ(8, tiger_it.Next<Tiger>()->ID());
    EXPECT_EQ(8, tiger_it.Get<Tiger>()->ID());
    EXPECT_EQ((Tiger*)(0), tiger_it.Next<Tiger>());
    EXPECT_EQ((Tiger*)(0), tiger_it.Get<Tiger>());

    /* Iterates through all Tiger objects, except id 5, which is deleted. */
    EXPECT_TRUE(idx.Delete(5));
    ObjIdx<Feline>::Iterator tiger_it_1 = idx.Begin();
    EXPECT_EQ(2, tiger_it_1.Next<Tiger>()->ID());
    EXPECT_EQ(2, tiger_it_1.Get<Tiger>()->ID());
    EXPECT_EQ(8, tiger_it_1.Next<Tiger>()->ID());
    EXPECT_EQ(8, tiger_it_1.Get<Tiger>()->ID());
    EXPECT_EQ((Tiger*)(0), tiger_it_1.Next<Tiger>());
    EXPECT_EQ((Tiger*)(0), tiger_it_1.Get<Tiger>());

    //idx.Print<Tiger>();
    //idx.Print();

    EXPECT_TRUE(idx.Delete(9));
    EXPECT_TRUE(idx.Delete(8));
    EXPECT_TRUE(idx.Delete(7));
    EXPECT_TRUE(idx.Delete(6));
    /* ID 5 was already deleted. */
    EXPECT_FALSE(idx.Delete(5));
    EXPECT_TRUE(idx.Delete(4));
    EXPECT_TRUE(idx.Delete(3));
    EXPECT_TRUE(idx.Delete(2));
    EXPECT_TRUE(idx.Delete(1));
    EXPECT_TRUE(idx.Delete(0));

    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);
  }

  TEST(ObjIdx, recycling_ids) {
    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);

    ObjIdx<Feline> idx;

    for (int i=0; i<100; i++){ EXPECT_EQ(i, idx.Create()->ID()); }
    EXPECT_EQ(100, idx.GetSize());
    EXPECT_EQ(100, idx.GetMaxCt());
    EXPECT_EQ(100, Feline::s_ct);

    for (int i=0; i<50; i++){ EXPECT_TRUE(idx.Delete(i)); }
    EXPECT_EQ(50, idx.GetSize());
    EXPECT_EQ(100, idx.GetMaxCt());
    EXPECT_EQ(50, Feline::s_ct);

    /* First 50 IDs should now be recycled. */
    for (int i=49; 0<=i; i--){ EXPECT_EQ(i, idx.Create()->ID()); }
    EXPECT_EQ(100, idx.GetSize());
    EXPECT_EQ(100, idx.GetMaxCt());
    EXPECT_EQ(100, Feline::s_ct);

    /* Next 50 IDs should start from 100. */
    for (int i=100; i<150; i++){ EXPECT_EQ(i, idx.Create()->ID()); }
    EXPECT_EQ(150, idx.GetSize());
    EXPECT_EQ(150, idx.GetMaxCt());
    EXPECT_EQ(150, Feline::s_ct);

    for (int i=0; i<150; i++){ EXPECT_TRUE(idx.Delete(i)); }
    EXPECT_EQ(0, idx.GetSize());
    EXPECT_EQ(150, idx.GetMaxCt());
    EXPECT_EQ(0, Feline::s_ct);

    EXPECT_EQ(0, Feline::s_ct);
    EXPECT_EQ(0, Lion::s_ct);
    EXPECT_EQ(0, Tiger::s_ct);
  }
}


//TEST(cStrand, brainstorm_workflow){
//  /* Load initial organism string */
////  Genome genome;
////  
////  cString config_filename = "avida-fsm.cfg";
////  cString organism_filename = "default-fsm.org";
////  cUserFeedbackReporter feedback_reporter;
////  tDictionary<cString> defs;
////  cAvidaConfig* cfg = new cAvidaConfig();
////  
////  // Initialize the configuration data...
////  chdir("work");
////  Avida::Initialize();
////  cfg->Load(config_filename, cString(Apto::FileSystem::GetCWD()), &feedback_reporter.fb(), &defs, /* flag_warn_default */ false);
////  cWorld* world = cWorld::Initialize(cfg, cString(Apto::FileSystem::GetCWD()), &feedback_reporter.fb());
////  feedback_reporter.ReportNewFeedback("Feedback from cWorld initialization:");
////  genome.LoadFromDetailFile(organism_filename, world->GetWorkingDir(), world->GetHardwareManager(), feedback_reporter.fb());
////  feedback_reporter.ReportNewFeedback("Feedback from genome.LoadFromDetailFile:");
////  chdir("..");
////  
////  cAvidaContext& ctx = world->GetDefaultContext();
////  cOrganism* organism = new cOrganism(world, ctx, genome, -1, SRC_ORGANISM_FILE_LOAD);
////  cHardwareFSM *hardware = dynamic_cast<cHardwareFSM*>(&organism->GetHardware());
////  EXPECT_EQ(HARDWARE_TYPE_CPU_FSM, hardware->GetType());
//  
//  /* I wish I had a... strand index to load the first string. */
//  /* Create an initial finite state machine. */
//  /* load initial string into index. */
//  // cStrand &strand(cStrandIndex.Add(genome.GetSequence()));
//
//  /*
//   I wish I had a library of strands; each strand references a string, and multiple strands can reference the same string. Each strand is inherently associated with binding sites, and possibly with FSMs.
//   I wish I had a library of binding sites; each binding sites has an ID, and multiple binding sites can have the same ID. Each label can be associated with a strand and position, or an FSM and a position.
//   I wish I had a library of FSMs; each FSM has program, and multiple FSMs can reference the same program; the program can be hard-coded, or can be encoded by a strand. Each FSM has a nonempty set of binding sites.
//   I wish I had a library of bindings; each binding references a set of binding criteria, and multiple bindings can reference the same criteria. Each binding consists of two or more binding sites.
//   */
//  
//  /*
//   Things that can happen:
//   - FSMs can execute
//   - Bindings can form
//   - Bindings can break
//   */
//  
//  /* I think that the bootstrap FSM depends on the contents of the index! */
//  /* I think the index can initially have multiple strands, if the initial sequence so encodes. */
//  /* I think that the bootstrap FSM can decode the initial sequence into multiple strands! */
//  /*
//   I think that the initial configuration can specify a collection of initial FSMs.
//   From each available FSMs, we get a list of binding sites.
//   */
//  /* make new state machine from strand... */
//  // cStrandFSM
//  
////  cStrandIndex strand_idx;
//  /* bind state machine to strand. */
//  
////  EXPECT_EQ(0, strand_idx.m_strands.GetSize());
////  strand_idx.Add(*strand, hardware);
////  EXPECT_EQ(1, strand_idx.m_strands.GetSize());
////  EXPECT_TRUE(strand_idx.Remove(*strand));
////  strand_idx.RemoveFromIndex(*strand);
////  EXPECT_EQ(0, strand_idx.m_strands.GetSize());
////  EXPECT_FALSE(strand_idx.Remove(*strand));
////
////  delete strand;
//
////  cStrand *strand = 0;
//
//}

//int cFSMDB::InsertSequence(const Apto::String &sequence) {
  //if (m_seqs.Has(seq)) { return m_seqs.Insert(seq); }

  ///* Scan new sequence for labels. */
  //int seq_id = m_seqs.Insert(seq);
  //Apto::Array<cHit, Apto::Smart> hits(bScanForLabels(seq_id, *this));
  //for (int i=0; i<hits.GetSize(); i++) {
  //  //LinkSeqLblPos(seq_id, hits[i].Lbl(), hits[i].Pos());
  //}
  //
  //return seq_id;
  //return -1;
//}
