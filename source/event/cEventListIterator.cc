//////////////////////////////////////////////////////////////////////////////
// Copyright (C) 1993 - 2003 California Institute of Technology             //
//                                                                          //
// Read the COPYING and README files, or contact 'avida@alife.org',         //
// before continuing.  SOME RESTRICTIONS MAY APPLY TO USE OF THIS FILE.     //
//////////////////////////////////////////////////////////////////////////////

#ifndef EVENT_LIST_ITERATOR_HH
#include "cEventListIterator.h"
#endif

#ifndef EVENT_LIST_HH
#include "cEventList.h"
#endif

using namespace std;

/////////////////
//  cEventListIterator
//  added by kaben.
/////////////////

void
cEventListIterator::PrintEvent(ostream& os){
  cEventList::PrintEvent(m_node, os);
}

