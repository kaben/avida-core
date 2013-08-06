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
#include "avida/core/InstructionSequence.h"
#include "apto/core/Array.h"
#include "apto/core/String.h"
#include "apto/core/StringBuffer.h"
#include "apto/core/Map.h"

class IDBase;
class ObjBase;
class cBasicLabelUtils;
class cFSM;
class cFSMDB;
class cFSMDef;
class cLabel;
class cLabelHit;
class cLabelIdx;
class cSeqIdx;
class cSequence;
class cStrand;
template <class T> class ObjIdx;

Apto::Array<cLabelHit, Apto::Smart> bScanForLabels(const Apto::String&, const cBasicLabelUtils&);


class IDBase {
protected:
  unsigned int m_obj_ct;
  unsigned int m_next_new_id;
  Apto::Array<int> m_recycled_ids;
protected:
  IDBase(): m_obj_ct(0), m_next_new_id(0) {}
  int NextID();
  void RecycleID(int id);
public:
  int GetSize() { return m_obj_ct; }
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


class cFSM : public ObjBase {
public:
  int m_fsm_def_id;
};


class cFSMDef : public ObjBase {
public:
  Apto::Map<int, Apto::Set<int> > m_lbl_sites;
  Apto::Set<int> m_fsm_ids;
};


class cLabel : public ObjBase {
public:
  Apto::Set<int> m_seq_ids;
  Apto::Set<int> m_fsm_def_ids;
};


class cLabelHit {
public:
  int m_label_id;
  int m_pos;
public:
  cLabelHit(int label_id = -1, int pos = -1): m_label_id(label_id), m_pos(pos) {}
public:
  int Lbl() const { return m_label_id; }
  int Pos() const { return m_pos; }
  bool operator==(const cLabelHit & other) const { return (other.m_label_id == m_label_id) && (other.m_pos == m_pos); }
  bool operator!=(const cLabelHit & other) const { return (other.m_label_id != m_label_id) || (other.m_pos != m_pos); }
};


class cLabelIdx {
protected:
  Apto::Map<int, cLabel> m_id2obj;
public:
  typedef typename Apto::Map<int, cLabel> AM;
  typedef typename AM::KeyIterator IDIter;
  typedef typename AM::ValueIterator ObjIter;
public:
  bool Has(int id);
  cLabel* Get(int id);
  int GetSize();
  cLabel* Insert(int id);
  bool Delete(int id);
  IDIter IDs();
  ObjIter Objs();
};


class cSeqIdx : public IDBase {
protected:
  Apto::Map<int, cSequence> m_id2obj;
  Apto::Map<Apto::String, int> m_str2id;
  Apto::Map<int, Apto::String> m_id2str;
public:
  typedef typename Apto::Map<int, cSequence> AM;
  typedef typename AM::KeyIterator IDIter;
  typedef typename AM::ValueIterator ObjIter;
public:
  bool Has(int id);
  bool Has(const Apto::String &str);
  cSequence* Get(int id);
  int GetID(const Apto::String &str);
  Apto::String GetString(int id);
  cSequence* Insert(const Apto::String &str);
  bool Delete(int id);
  bool Delete(const Apto::String &str);
  IDIter IDs();
  ObjIter Objs();
};


class cSequence : public ObjBase {
public:
  Apto::Map<int, Apto::Set<int> > m_lbl_sites;
  Apto::Set<int> m_strand_ids;
};


class cStrand : public ObjBase {
public:
  int m_seq_id;
};


template <class T>
class ObjIdx : public IDBase {
protected:
  Apto::Map<int, T> m_id2obj;
public:
  typedef typename Apto::Map<int, T> AM;
  typedef typename AM::KeyIterator IDIter;
  typedef typename AM::ValueIterator ObjIter;
public:
  bool Has(int id) { return m_id2obj.Has(id); }
  T* Get(int id) { return (Has(id))?(&m_id2obj.Get(id)):(0); }
  T* Create() {
    int id = NextID();
    T* ptr = &m_id2obj.Get(id);
    ptr->m_id = id;
    return ptr;
  }
  bool Delete(int id) {
    if (Has(id)) {
      m_id2obj.Remove(id);
      RecycleID(id);
      return true;
    }
    return false;
  }
  IDIter IDs() { return m_id2obj.Keys(); }
  ObjIter Objs() { return m_id2obj.Values(); }
};


class cFSMDB {
public:
  cLabelIdx m_lbls;
  cSeqIdx m_seqs;
  ObjIdx<cStrand> m_strands;
  ObjIdx<cFSMDef> m_fsm_defs;
  ObjIdx<cFSM> m_fsms;

  cBasicLabelUtils m_label_utils;

  int CreateStrand(const Apto::String &sequence);
  bool RemoveStrand(int strand_id);
protected:
  void LinkSeqLblPos(int seq_id, int lbl_id, int lbl_pos);
  int InsertSequence(const Apto::String &sequence);
  void UnlinkSeqLbls(int seq_id);
  bool RemoveSequence(int seq_id);
};


#endif
