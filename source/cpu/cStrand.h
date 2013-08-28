/*
 *  cStrand.h
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

#ifndef cStrand_h
#define cStrand_h

#include "cCodeLabel.h"
#include "cCPUMemory.h"
#include "avida/core/InstructionSequence.h"
#include "apto/core/Array.h"
#include "apto/core/Functor.h"
#include "apto/core/Map.h"
#include "apto/core/PriorityScheduler.h"
#include "apto/core/PriorityScheduler.h"
#include "apto/core/Random.h"
#include "apto/core/SmartPtr.h"
#include "apto/core/String.h"
#include "apto/core/StringBuffer.h"
#include "apto/rng/AvidaRNG.h"

#include <iostream>


class IDs;
class ObjBase;
class cBasicLabelUtils;
class cBindable;
class cFSM;
class cFSMFunctorObject;
class cFSMDB;
class cFSMDef;
class cHit;
class cLabel;
class cLabelDeletemeIdx;
class cSeqIdx;
class cSequence;
class cStrand;
template <class T> class ObjIdx;

Apto::Array<cHit, Apto::Smart> bScanForLabels(const Apto::String&, const cBasicLabelUtils&);


typedef Apto::Functor<void, Apto::TL::Create<int> > FSMFunctor;


namespace Apto {
  namespace Scheduler {
    class ProbabilisticDynamic : public PriorityScheduler {
    protected:
      SmartPtr<Random> m_rng;
    public:
      inline ProbabilisticDynamic(int num_entries, SmartPtr<Random> rng) : m_rng(rng), m_index(num_entries) { ; }
      ~ProbabilisticDynamic();
      void AdjustPriority(int entry_id, double priority);
      int Next();
      int GetSize() { return m_index.GetSize(); }
      void Resize(int size) { m_index.Resize(size); }
    protected:
      class DynamicWeightedIndex {
      public:
        int m_size;
        Array<double, Smart> m_item_weight;
        Array<double, Smart> m_subtree_weight;
      public:
        DynamicWeightedIndex(int size);
        ~DynamicWeightedIndex() { ; }
        void SetWeight(int entry_id, double weight);
        double TotalWeight() const { return m_subtree_weight[0]; }
        int FindPosition(double position) { return findPosition(position, 0); }
        int GetSize() { return m_size; }
        void Resize(int size);
      public:
        int findPosition(double position, int root_id);
        int parentOf(int entry_id) { return (entry_id - 1) / 2; }
        int leftChildOf(int entry_id) { return 2 * entry_id + 1; }
        int rightChildOf(int entry_id) { return 2 * entry_id + 2; }
      };
      DynamicWeightedIndex m_index;
    };

    class ProbabilisticID : public ProbabilisticDynamic {
    public:
      Map<int, int> id2pos;
      Array<int, Smart> pos2id;
      inline ProbabilisticID(int num_entries, SmartPtr<Random> rng)
      : ProbabilisticDynamic(num_entries, rng)
      , pos2id(num_entries)
      {}
      void AdjustIDPriority(int id, double priority);
      int NextID();
    };

    // Integrated
    // --------------------------------------------------------------------------------------------------------------
    //  The integrated scheduler method relies on breaking up all merits into sums of powers of 2 (i.e. using the
    //  binary representation of the priority). All items with merits in the highest power of two will get the most
    //  time, and subsequent priority components will have time divided, continuing recursively. The simplest way of
    //  doing this while maximizing evenness of distribution of time slices is to simply alternate executing the best,
    //  and everything else (where in everything else we again alternate with the best of this sub-list recursively).
    class IntegratedDynamic : public PriorityScheduler {
    public:
      struct Node;
    public:
      Array<Node*> m_node_array;
      Array<double, Smart> m_priority_chart;
    private:
      IntegratedDynamic();
      IntegratedDynamic(const IntegratedDynamic&);
      IntegratedDynamic& operator=(const IntegratedDynamic&);
    public:
      IntegratedDynamic(int entry_count);
      ~IntegratedDynamic();
      void AdjustPriority(int entry_id, double priority);
      int Next();
      int GetSize() { return m_priority_chart.GetSize(); }
      void Resize(int size);
    public:
      void insertNode(int node_id);
      void removeNode(int node_id);
      void resizeNodes(int new_size);
    public:
      struct Node {
        /*
        Each entry in this array corresponds to the item with the same ID. if
        the entry is not in the list, its value in the array will be 0. If it
        is in the list, it will point to the cell of the next included entry.
        The last included entry has a -1 in its cell.
        */
        Array<int, Smart> active_array; 
        int first_entry;          // ID of the first active entry
        int active_entry;         // ID of the next scheduled entry
        int node_id;              // A unique ID (representative the relative priority bit
        
        int size;                 // Number of active items in this node
        int process_size;         // Number of times this node should be executed before the next node is.
        int process_count;        // Number of times this node has been executed
        bool execute;             // Should this node execute or pass?
        
        Node* next;
        Node* prev;
        
        inline Node(int entry_count = 0, int in_node_id = -1)
        : active_array(entry_count)
        , first_entry(-1)
        , active_entry(-1)
        , node_id(in_node_id)
        , size(0)
        , process_size(1)
        , process_count(0)
        , execute(true)
        , next(NULL)
        , prev(NULL)
        { active_array.SetAll(0); }
        inline ~Node() {}
        void Insert(int entry_id);
        void Remove(int entry_id);
        int Next();
      };
    };

    class IntegratedID : public IntegratedDynamic {
    public:
      Map<int, int> id2pos;
      Array<int, Smart> pos2id;
      inline IntegratedID(int num_entries)
      : IntegratedDynamic(num_entries)
      , pos2id(num_entries)
      {}
      void AdjustIDPriority(int id, double priority);
      int NextID();
    };
  }
}


class IDMgr {
protected:
  unsigned int m_obj_ct;
  unsigned int m_max_ct;
  Apto::Array<int> m_recycled_ids;
public:
  IDMgr(): m_obj_ct(0), m_max_ct(0) {}
  int NextID();
  void RecycleID(int id);
public:
  int GetSize() { return m_obj_ct; }
  int GetMaxCt() { return m_max_ct; }
};


class ObjBase {
  template <class T> friend class ObjIdx;
  friend class cSeqIdx;
  friend class cLabelIdx;
private:
  int m_id;
public:
  int ID() const { return m_id; }
};


class cBasicLabelUtils {
public:
  char m_zero_char;
  int m_num_nops;
  int m_max_label_size;
public:
  int GetMaxLabelSize() const { return m_max_label_size; }
  int GetNumNops() const { return m_num_nops; }
  int IsNop(char c) const;
  int Seq2ID(const Apto::String& seq) const;
  Apto::StringBuffer ID2Seq(const int id) const;
  int SeqLen(const int id) const;
  int Rotate(const int id, const int rot = 2) const;
  int ReverseComplement(const int id, const int rot = 2) const;
public:
  cBasicLabelUtils(
    const char zero_char = 'a',
    const int num_nops = 4,
    const int max_label_size = 3
  )
  : m_zero_char(zero_char)
  , m_num_nops(num_nops)
  , m_max_label_size(max_label_size)
  {}
};


class cHalfBinding : public ObjBase {
public:
  int m_parent_id;
  int m_lbl_id;
  int m_lbl_pos;
  int m_lbl_len;
  int m_other_half_binding_id;
  void Set(int parent_id, int lbl_id, int lbl_pos, int lbl_len, int other_half_binding_id){
    m_parent_id = parent_id;
    m_lbl_id = lbl_id;
    m_lbl_pos = lbl_pos;
    m_lbl_len = lbl_len;
    m_other_half_binding_id = other_half_binding_id;
  }
};


class cHit {
public:
  int m_id;
  int m_pos;
public:
  cHit(int id = -1, int pos = -1): m_id(id), m_pos(pos) {}
public:
  int Lbl() const { return m_id; }
  int Pos() const { return m_pos; }
  bool operator==(const cHit & other) const { return (other.m_id == m_id) && (other.m_pos == m_pos); }
  bool operator!=(const cHit & other) const { return (other.m_id != m_id) || (other.m_pos != m_pos); }
};

class cLabel : public ObjBase {
public:
  Apto::Set<int> m_seq_ids;
  Apto::Set<int> m_fsm_def_ids;
};

class cLabeled : public ObjBase {
public:
  /* Index by label ID of label positions. */
  Apto::Map<int, Apto::Array<int> > m_labels;
public:
  virtual ~cLabeled(){}
};

class cSequence : public cLabeled {
public:
  Apto::Set<int> m_strand_ids;
};

class cFSMDef : public cLabeled {
public:
  Apto::Set<int> m_fsm_ids;
public:
  virtual ~cFSMDef(){}
};

class cFunction {
public:
};

class cFSMBootstrapDef : public cFSMDef {
public:
};

class cNFADef : public cFSMDef {
public:
  /* Use: int next_state_id = m_transition_relation[state_id][symbol_id][next_state_idx]; */
  Apto::Map<int, Apto::Map<int, Apto::Array<int> > > m_transition_relation;
  /* Use: int function_id = m_function_relation[state_id][function_idx]; */
  Apto::Map<int, Apto::Array<int> > m_function_relation;
};

class cBindable: public ObjBase {
public:
  /* Index by position of bindsites. */
  Apto::Map<int, Apto::Set<int> > m_bindpts;
public:
  virtual ~cBindable(){}
  Apto::Array<int> GetBindings(cFSMDB &db);
  virtual Apto::Map<int, Apto::Array<int> > &GetLabels(cFSMDB &db) = 0;
  virtual Apto::String AsString(cFSMDB &db) = 0;
};

class cStrand : public cBindable {
public:
  int m_seq_id;
  virtual Apto::Map<int, Apto::Array<int> > &GetLabels(cFSMDB &db);
  virtual Apto::String AsString(cFSMDB &db);
};

class cFSM : public cBindable {
public:
  int m_fsm_def_id;
  virtual Apto::Map<int, Apto::Array<int> > &GetLabels(cFSMDB &db);
  virtual Apto::String AsString(cFSMDB &db);
};

class cFSMBootstrap : public cFSM {
public:
};

class cNFA : public cFSM {
public:
  int m_state_id;
  Apto::SmartPtr<Apto::RNG::AvidaRNG> m_rng;
public:
  int Transition(int symbol_id, cFSMDB &db);
};

template <class T>
class ObjIdx : public IDMgr {
protected:
  Apto::Array<T*, Apto::Smart> m_objs;
public:
  class Iterator;
  friend class Iterator;
public:
  ~ObjIdx() { for (Iterator it = Begin(); it.Next();) { Delete(it.ID()); } }
public:
  /*
  This is a polymorphic factory function using the default constructor for U.
  It will only work if U and T are the same type, or if U is a polymorphic
  subclass of T.

  Requires U to have default constructor and assignment operator.
  */
  template <class U>
  U* Create() {
    int id = NextID();
    if (m_objs.GetSize() < GetMaxCt()) { m_objs.Resize(GetMaxCt()); }
    U* ptr = new U;
    ptr->m_id = id;
    m_objs[id] = ptr;
    return ptr;
  }
  /*
  Requires T to have default constructor and assignment operator.
  */
  T* Create() { return Create<T>(); }
  /*
  If id < GetMaxCt() then either id is a valid index into m_objs, or it was
  once a valid index into m_objs. If it was once but no longer is a valid
  index, then m_objs[id] was deleted and should now return a null pointer. If
  it is still a valid index, then m_objs[id] should return a non-null pointer
  to a valid object.
  */
  bool Has(int id) {
    return (0 <= id) && (id < GetMaxCt()) && (m_objs[id] != (T*)(0));
  }
  /*
  This is a polymorphic getter. If the indexed object cannot be dynamically
  cast to type U, then the getter will return a null pointer.
  */
  template <class U> U* Get(int id) {
    return ((0 <= id) && (id < GetMaxCt()))?(dynamic_cast<U*>(m_objs[id])):((U*)(0));
  }
  T* Get(int id) {
    return ((0 <= id) && (id < GetMaxCt()))?(m_objs[id]):((T*)(0));
  }
  bool Delete(int id) {
    if (Has(id)) {
      delete m_objs[id];
      m_objs[id] = (T*)(0);
      RecycleID(id);
      return true;
    }
    return false;
  }
  Iterator Begin() { return Iterator(*this); }
  template <class U> void Print() {
    for (Iterator it = Begin(); it.Next();) {
      U* ptr = it.Get<U>();
      if (ptr) { std::cout << "id: " << it.ID() << ", ptr: " << ptr << std::endl; }
    }
  }
  void Print() { Print<T>(); }
public:
  class Iterator {
    friend class ObjIdx;
  protected:
    ObjIdx& m_idx;
    int m_id;
    Iterator(); // @not_implemented
    Iterator(ObjIdx& idx) : m_idx(idx), m_id(-1) { ; }
  public:      
    int ID() { return m_id; }
    template <class U>
    U* Get() { return m_idx.Get<U>(m_id); }
    T* Get() { return m_idx.Get(m_id); }
    template <class U>
    U* Next() {
      while (++m_id < m_idx.GetMaxCt()) {
        U* ptr = dynamic_cast<U*>(m_idx.m_objs[m_id]);
        if (ptr != (U*)(0)) { return ptr; }
      }
      return (U*)(0);
    }
    T* Next() {
      while (++m_id < m_idx.GetMaxCt()) {
        T* ptr = m_idx.m_objs[m_id];
        if (ptr != (T*)(0)) { return ptr; }
      }
      return (T*)(0);
    }
  };
};


class cSeqIdx : public IDMgr {
protected:
  Apto::Array<cSequence*, Apto::Smart> m_objs;
  Apto::Map<Apto::String, int> m_str2id;
  Apto::Map<int, Apto::String> m_id2str;
public:
  ~cSeqIdx() {
    /* FIXME: cleanup. */
    //Apto::Array<int, Apto::Smart> ids;
    //for (Apto::Map<int, Apto::String>::KeyIterator it = m_id2str.Keys(); it.Next();) { ids.Push(*it.Get()); }
    //for (int i=0; i<ids.GetSize(); i++){ Delete(ids[i]); }
  }
public:
  bool Has(int id);
  bool Has(const Apto::String &str);
  cSequence* Get(int id);
  int GetID(const Apto::String &str);
  Apto::String GetString(int id);
  cSequence* Insert(const Apto::String &str);
  bool Delete(int id);
  bool Delete(const Apto::String &str);
};


class cLabelIdx : public IDMgr {
public:
  class Iterator;
  friend class Iterator;
protected:
  Apto::Array<cLabel*, Apto::Smart> m_objs;
  Apto::Map<int, int> m_id2idx;
public:
  ~cLabelIdx() {
    /* FIXME: cleanup. */
    //Apto::Array<int, Apto::Smart> ids;
    //for (Iterator it = Begin(); it.Next();) { ids.Push(it.ID()); }
    //for (int i=0; i<ids.GetSize(); i++){ Delete(ids[i]); }
  }
public:
  bool Has(int id);
  cLabel* Get(int id);
  cLabel* Insert(int id);
  bool Delete(int id);
  Iterator Begin() { return Iterator(*this); }
public:
  class Iterator {
    friend class cLabelIdx;
  protected:
    cLabelIdx& m_idx;
    Apto::Map<int, int>::KeyIterator m_it;
    Iterator(); // @not_implemented
    Iterator(cLabelIdx& idx);
  public:      
    int ID();
    cLabel* Get();
    cLabel* Next();
  };
};


class cKinetics {
public:
  Apto::SmartPtr<Apto::RNG::AvidaRNG> m_rng;
  Apto::Scheduler::ProbabilisticDynamic m_collision_scheduler;
  ObjIdx<cBindable> m_bindables;
  ObjIdx<cHalfBinding> m_half_bindings;
public:
  cKinetics(Apto::SmartPtr<Apto::RNG::AvidaRNG> rng)
  : m_rng(rng)
  , m_collision_scheduler(0, m_rng)
  {};
};

class cFSMFunctorObject {
  cFSMDB &m_db;
public:
  cFSMFunctorObject(cFSMDB &db);
public:
  void Function0(int caller_id);
  void Function1(int caller_id);
  void Function2(int caller_id);
  void Function3(int caller_id);
  void Function4(int caller_id);
  void Function5(int caller_id);
  void Function6(int caller_id);
  void Function7(int caller_id);
  void Function8(int caller_id);
  void Function9(int caller_id);
};

class cFSMDB {
public:
  Apto::SmartPtr<Apto::RNG::AvidaRNG> m_rng;
  cBasicLabelUtils m_label_utils;

  cLabelIdx m_lbls;
  cSeqIdx m_seqs;
  ObjIdx<cFSMDef> m_fsm_defs;
  ObjIdx<cBindable> m_bindables;
  ObjIdx<cHalfBinding> m_half_bindings;
  //cKinetics m_kinetics;
  Apto::Scheduler::ProbabilisticDynamic m_collision_scheduler;
  Apto::Scheduler::IntegratedDynamic m_unbinding_scheduler;

  Apto::Array<FSMFunctor> m_functors;
  cFSMFunctorObject m_functor_object;

  int CreateStrand();
  int CreateStrand(const Apto::String &sequence);
  void AssociateSeqToStrand(int strand_id, const Apto::String &seq);
  void RemoveStrand(int strand_id);
  void LyseStrand(int strand_id, int at_position);

  int CreateFSMBootstrap();

  bool SingleCollision();
  bool Collide(int bindable_id_0, int bindable_id_1);
  bool SingleUnbinding();
  bool Unbind(int half_binding_id);
  bool SingleRebinding();

public:
  cFSMDB();
protected:
  int InsertSequence(const Apto::String &sequence);
  void UnlinkSeqLbls(int seq_id);
  bool RemoveSequence(int seq_id);
};


#endif
