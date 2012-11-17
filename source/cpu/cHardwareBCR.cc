/*
 *  cHardwareBCR.cc
 *  Avida
 *
 *  Created by David on 11/2/2012 based on cHardwareMBE.h
 *  Copyright 2012 Michigan State University. All rights reserved.
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
 *  Authors: David M. Bryson <david@programerror.com>, Aaron P. Wagner <apwagner@msu.edu>
 *
 */


#include "cHardwareBCR.h"

#include "avida/core/Feedback.h"
#include "avida/core/WorldDriver.h"

#include "cAvidaContext.h"
#include "cHardwareManager.h"
#include "cHardwareTracer.h"
#include "cInstSet.h"
#include "cOrganism.h"
#include "cPhenotype.h"
#include "cPopulation.h"
#include "cStateGrid.h"
#include "cWorld.h"

#include "tInstLibEntry.h"

#include <climits>
#include <fstream>

using namespace std;
using namespace Avida;
using namespace AvidaTools;



tInstLib<cHardwareBCR::tMethod>* cHardwareBCR::s_inst_slib = cHardwareBCR::initInstLib();

tInstLib<cHardwareBCR::tMethod>* cHardwareBCR::initInstLib(void)
{
  struct cNOPEntry {
    cString name;
    int nop_mod;
    cNOPEntry(const cString &name, int nop_mod)
    : name(name), nop_mod(nop_mod) {}
  };
  static const cNOPEntry s_n_array[] = {
    cNOPEntry("nop-A", rAX),
    cNOPEntry("nop-B", rBX),
    cNOPEntry("nop-C", rCX),
    cNOPEntry("nop-D", rDX),
    cNOPEntry("nop-E", rEX),
    cNOPEntry("nop-F", rFX),
    cNOPEntry("nop-G", rGX),
    cNOPEntry("nop-H", rHX),
    cNOPEntry("nop-I", rIX),
    cNOPEntry("nop-J", rJX),
    cNOPEntry("nop-K", rKX),
    cNOPEntry("nop-L", rLX),
    cNOPEntry("nop-M", rMX),
    cNOPEntry("nop-N", rNX),
    cNOPEntry("nop-O", rOX),
    cNOPEntry("nop-P", rPX),
  };
  
  static const tInstLibEntry<tMethod> s_f_array[] = {
    /*
     Note: all entries of cNOPEntryCPU s_n_array must have corresponding
     in the same order in tInstLibEntry<tMethod> s_f_array, and these entries must
     be the first elements of s_f_array.
     */
    tInstLibEntry<tMethod>("nop-A", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-B", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-C", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-D", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-E", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-F", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-G", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-H", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-I", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-J", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-K", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-L", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-M", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-N", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-O", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    tInstLibEntry<tMethod>("nop-P", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, nInstFlag::NOP, "No-operation; modifies other instructions"),
    
    tInstLibEntry<tMethod>("NULL", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, 0, "True no-operation instruction: does nothing"),
    tInstLibEntry<tMethod>("nop-X", &cHardwareBCR::Inst_Nop, INST_CLASS_NOP, 0, "True no-operation instruction: does nothing"),
    
    // Threading 
    tInstLibEntry<tMethod>("thread-create", &cHardwareBCR::Inst_ThreadCreate, INST_CLASS_OTHER, 0, "", BEHAV_CLASS_NONE),
    tInstLibEntry<tMethod>("thread-cancel", &cHardwareBCR::Inst_ThreadCancel, INST_CLASS_OTHER, 0, "", BEHAV_CLASS_NONE),
    tInstLibEntry<tMethod>("thread-id", &cHardwareBCR::Inst_ThreadID, INST_CLASS_OTHER, 0, "", BEHAV_CLASS_NONE),
    tInstLibEntry<tMethod>("yield", &cHardwareBCR::Inst_Yield, INST_CLASS_OTHER, 0, "", BEHAV_CLASS_NONE),

    // Standard Conditionals
    tInstLibEntry<tMethod>("if-n-equ", &cHardwareBCR::Inst_IfNEqu, INST_CLASS_CONDITIONAL, 0, "Execute next instruction if ?BX?!=?CX?, else skip it"),
    tInstLibEntry<tMethod>("if-less", &cHardwareBCR::Inst_IfLess, INST_CLASS_CONDITIONAL, 0, "Execute next instruction if ?BX? < ?CX?, else skip it"),
    tInstLibEntry<tMethod>("if-not-0", &cHardwareBCR::Inst_IfNotZero, INST_CLASS_CONDITIONAL, 0, "Execute next instruction if ?BX? != 0, else skip it"),
    tInstLibEntry<tMethod>("if-equ-0", &cHardwareBCR::Inst_IfEqualZero, INST_CLASS_CONDITIONAL, 0, "Execute next instruction if ?BX? == 0, else skip it"),
    tInstLibEntry<tMethod>("if-gtr-0", &cHardwareBCR::Inst_IfGreaterThanZero, INST_CLASS_CONDITIONAL, 0, "Execute next instruction if ?BX? > 0, else skip it"),
    tInstLibEntry<tMethod>("if-less-0", &cHardwareBCR::Inst_IfLessThanZero, INST_CLASS_CONDITIONAL, 0, "Execute next instruction if ?BX? < 0, else skip it"),
    tInstLibEntry<tMethod>("if-gtr-x", &cHardwareBCR::Inst_IfGtrX, INST_CLASS_CONDITIONAL),
    tInstLibEntry<tMethod>("if-equ-x", &cHardwareBCR::Inst_IfEquX, INST_CLASS_CONDITIONAL),
    
    // Core ALU Operations
    tInstLibEntry<tMethod>("pop", &cHardwareBCR::Inst_Pop, INST_CLASS_DATA, 0, "Remove top number from stack and place into ?BX?"),
    tInstLibEntry<tMethod>("push", &cHardwareBCR::Inst_Push, INST_CLASS_DATA, 0, "Copy number from ?BX? and place it into the stack"),
    tInstLibEntry<tMethod>("pop-all", &cHardwareBCR::Inst_PopAll, INST_CLASS_DATA, 0, "Remove top numbers from stack and place into ?BX?"),
    tInstLibEntry<tMethod>("push-all", &cHardwareBCR::Inst_PushAll, INST_CLASS_DATA, 0, "Copy number from all registers and place into the stack"),
    tInstLibEntry<tMethod>("swap-stk", &cHardwareBCR::Inst_SwitchStack, INST_CLASS_DATA, 0, "Toggle which stack is currently being used"),
    tInstLibEntry<tMethod>("swap-stk-top", &cHardwareBCR::Inst_SwapStackTop, INST_CLASS_DATA, 0, "Swap the values at the top of both stacks"),
    tInstLibEntry<tMethod>("swap", &cHardwareBCR::Inst_Swap, INST_CLASS_DATA, 0, "Swap the contents of ?BX? with ?CX?"),
    
    tInstLibEntry<tMethod>("shift-r", &cHardwareBCR::Inst_ShiftR, INST_CLASS_ARITHMETIC_LOGIC, 0, "Shift bits in ?BX? right by one (divide by two)"),
    tInstLibEntry<tMethod>("shift-l", &cHardwareBCR::Inst_ShiftL, INST_CLASS_ARITHMETIC_LOGIC, 0, "Shift bits in ?BX? left by one (multiply by two)"),
    tInstLibEntry<tMethod>("inc", &cHardwareBCR::Inst_Inc, INST_CLASS_ARITHMETIC_LOGIC, 0, "Increment ?BX? by one"),
    tInstLibEntry<tMethod>("dec", &cHardwareBCR::Inst_Dec, INST_CLASS_ARITHMETIC_LOGIC, 0, "Decrement ?BX? by one"),
    tInstLibEntry<tMethod>("zero", &cHardwareBCR::Inst_Zero, INST_CLASS_ARITHMETIC_LOGIC, 0, "Set ?BX? to 0"),
    tInstLibEntry<tMethod>("one", &cHardwareBCR::Inst_One, INST_CLASS_ARITHMETIC_LOGIC, 0, "Set ?BX? to 0"),
    tInstLibEntry<tMethod>("rand", &cHardwareBCR::Inst_Rand, INST_CLASS_ARITHMETIC_LOGIC, 0, "Set ?BX? to rand number"),
    
    tInstLibEntry<tMethod>("add", &cHardwareBCR::Inst_Add, INST_CLASS_ARITHMETIC_LOGIC, 0, "Add BX to CX and place the result in ?BX?"),
    tInstLibEntry<tMethod>("sub", &cHardwareBCR::Inst_Sub, INST_CLASS_ARITHMETIC_LOGIC, 0, "Subtract CX from BX and place the result in ?BX?"),
    tInstLibEntry<tMethod>("nand", &cHardwareBCR::Inst_Nand, INST_CLASS_ARITHMETIC_LOGIC, 0, "Nand BX by CX and place the result in ?BX?"),
    
    tInstLibEntry<tMethod>("IO", &cHardwareBCR::Inst_TaskIO, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "Output ?BX?, and input new number back into ?BX?", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("input", &cHardwareBCR::Inst_TaskInput, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "Input new number into ?BX?", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("output", &cHardwareBCR::Inst_TaskOutput, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "Output ?BX?", BEHAV_CLASS_ACTION),
    
    tInstLibEntry<tMethod>("mult", &cHardwareBCR::Inst_Mult, INST_CLASS_ARITHMETIC_LOGIC, 0, "Multiple BX by CX and place the result in ?BX?"),
    tInstLibEntry<tMethod>("div", &cHardwareBCR::Inst_Div, INST_CLASS_ARITHMETIC_LOGIC, 0, "Divide BX by CX and place the result in ?BX?"),
    tInstLibEntry<tMethod>("mod", &cHardwareBCR::Inst_Mod, INST_CLASS_ARITHMETIC_LOGIC),
        
    // Flow Control Instructions
    tInstLibEntry<tMethod>("label", &cHardwareBCR::Inst_Label, INST_CLASS_FLOW_CONTROL, nInstFlag::LABEL),
    tInstLibEntry<tMethod>("search-lbl-comp-s", &cHardwareBCR::Inst_Search_Label_Comp_S, INST_CLASS_FLOW_CONTROL, 0, "Find complement label from genome start and move the flow head"),
    tInstLibEntry<tMethod>("search-lbl-comp-f", &cHardwareBCR::Inst_Search_Label_Comp_F, INST_CLASS_FLOW_CONTROL, 0, "Find complement label forward and move the flow head"),
    tInstLibEntry<tMethod>("search-lbl-comp-b", &cHardwareBCR::Inst_Search_Label_Comp_B, INST_CLASS_FLOW_CONTROL, 0, "Find complement label backward and move the flow head"),
    tInstLibEntry<tMethod>("search-lbl-direct-s", &cHardwareBCR::Inst_Search_Label_Direct_S, INST_CLASS_FLOW_CONTROL, 0, "Find direct label from genome start and move the flow head"),
    tInstLibEntry<tMethod>("search-lbl-direct-f", &cHardwareBCR::Inst_Search_Label_Direct_F, INST_CLASS_FLOW_CONTROL, 0, "Find direct label forward and move the flow head"),
    tInstLibEntry<tMethod>("search-lbl-direct-b", &cHardwareBCR::Inst_Search_Label_Direct_B, INST_CLASS_FLOW_CONTROL, 0, "Find direct label backward and move the flow head"),
    tInstLibEntry<tMethod>("search-seq-comp-s", &cHardwareBCR::Inst_Search_Seq_Comp_S, INST_CLASS_FLOW_CONTROL, 0, "Find complement template from genome start and move the flow head"),
    tInstLibEntry<tMethod>("search-seq-comp-f", &cHardwareBCR::Inst_Search_Seq_Comp_F, INST_CLASS_FLOW_CONTROL, 0, "Find complement template forward and move the flow head"),
    tInstLibEntry<tMethod>("search-seq-comp-b", &cHardwareBCR::Inst_Search_Seq_Comp_B, INST_CLASS_FLOW_CONTROL, 0, "Find complement template backward and move the flow head"),
    tInstLibEntry<tMethod>("search-seq-direct-s", &cHardwareBCR::Inst_Search_Seq_Direct_S, INST_CLASS_FLOW_CONTROL, 0, "Find direct template from genome start and move the flow head"),
    tInstLibEntry<tMethod>("search-seq-direct-f", &cHardwareBCR::Inst_Search_Seq_Direct_F, INST_CLASS_FLOW_CONTROL, 0, "Find direct template forward and move the flow head"),
    tInstLibEntry<tMethod>("search-seq-direct-b", &cHardwareBCR::Inst_Search_Seq_Direct_B, INST_CLASS_FLOW_CONTROL, 0, "Find direct template backward and move the flow head"),

    tInstLibEntry<tMethod>("mov-head", &cHardwareBCR::Inst_MoveHead, INST_CLASS_FLOW_CONTROL, 0, "Move head ?IP? to the flow head"),
    tInstLibEntry<tMethod>("mov-head-if-n-equ", &cHardwareBCR::Inst_MoveHeadIfNEqu, INST_CLASS_FLOW_CONTROL, 0, "Move head ?IP? to the flow head if ?BX? != ?CX?"),
    tInstLibEntry<tMethod>("mov-head-if-less", &cHardwareBCR::Inst_MoveHeadIfLess, INST_CLASS_FLOW_CONTROL, 0, "Move head ?IP? to the flow head if ?BX? != ?CX?"),
    
    tInstLibEntry<tMethod>("jmp-head", &cHardwareBCR::Inst_JumpHead, INST_CLASS_FLOW_CONTROL, 0, "Move head ?Flow? by amount in ?CX? register"),
    tInstLibEntry<tMethod>("get-head", &cHardwareBCR::Inst_GetHead, INST_CLASS_FLOW_CONTROL, 0, "Copy the position of the ?IP? head into ?CX?"),

    tInstLibEntry<tMethod>("set-memory", &cHardwareBCR::Inst_SetMemory, INST_CLASS_FLOW_CONTROL, 0, "Set ?mem_space_label? of the ?Flow? head."),
    
    // Replication Instructions
    tInstLibEntry<tMethod>("divide", &cHardwareBCR::Inst_Divide, INST_CLASS_LIFECYCLE, nInstFlag::STALL, "Divide code between read and write heads.", BEHAV_CLASS_COPY),
    tInstLibEntry<tMethod>("divide-memory", &cHardwareBCR::Inst_DivideMemory, INST_CLASS_LIFECYCLE, nInstFlag::STALL, "Divide memory space.", BEHAV_CLASS_COPY),
    tInstLibEntry<tMethod>("h-copy", &cHardwareBCR::Inst_HeadCopy, INST_CLASS_LIFECYCLE, 0, "Copy from read-head to write-head; advance both", BEHAV_CLASS_COPY),
    tInstLibEntry<tMethod>("h-read", &cHardwareBCR::Inst_HeadRead, INST_CLASS_LIFECYCLE, 0, "Read instruction from ?read-head? to ?AX?; advance the head.", BEHAV_CLASS_COPY),
    tInstLibEntry<tMethod>("h-write", &cHardwareBCR::Inst_HeadWrite, INST_CLASS_LIFECYCLE, 0, "Write to ?write-head? instruction from ?AX?; advance the head.", BEHAV_CLASS_COPY),
    tInstLibEntry<tMethod>("if-copied-lbl-comp", &cHardwareBCR::Inst_IfCopiedCompLabel, INST_CLASS_CONDITIONAL, 0, "Execute next if we copied complement of attached label"),
    tInstLibEntry<tMethod>("if-copied-lbl-direct", &cHardwareBCR::Inst_IfCopiedDirectLabel, INST_CLASS_CONDITIONAL, 0, "Execute next if we copied direct match of the attached label"),
    tInstLibEntry<tMethod>("if-copied-seq-comp", &cHardwareBCR::Inst_IfCopiedCompSeq, INST_CLASS_CONDITIONAL, 0, "Execute next if we copied complement of attached sequence"),
    tInstLibEntry<tMethod>("if-copied-seq-direct", &cHardwareBCR::Inst_IfCopiedDirectSeq, INST_CLASS_CONDITIONAL, 0, "Execute next if we copied direct match of the attached sequence"),
    tInstLibEntry<tMethod>("did-copy-lbl-comp", &cHardwareBCR::Inst_DidCopyCompLabel, INST_CLASS_OTHER, 0, "Execute next if we copied complement of attached label"),
    tInstLibEntry<tMethod>("did-copy-lbl-direct", &cHardwareBCR::Inst_DidCopyDirectLabel, INST_CLASS_OTHER, 0, "Execute next if we copied direct match of the attached label"),
    tInstLibEntry<tMethod>("did-copy-seq-comp", &cHardwareBCR::Inst_DidCopyCompSeq, INST_CLASS_OTHER, 0, "Execute next if we copied complement of attached sequence"),
    tInstLibEntry<tMethod>("did-copy-seq-direct", &cHardwareBCR::Inst_DidCopyDirectSeq, INST_CLASS_OTHER, 0, "Execute next if we copied direct match of the attached sequence"),
    
    tInstLibEntry<tMethod>("repro", &cHardwareBCR::Inst_Repro, INST_CLASS_LIFECYCLE, nInstFlag::STALL, "Instantly reproduces the organism", BEHAV_CLASS_COPY),
    
    tInstLibEntry<tMethod>("die", &cHardwareBCR::Inst_Die, INST_CLASS_LIFECYCLE, nInstFlag::STALL, "Instantly kills the organism", BEHAV_CLASS_COPY),
    
    // Thread Execution Control
    tInstLibEntry<tMethod>("wait-cond-equ", &cHardwareBCR::Inst_WaitCondition_Equal, INST_CLASS_OTHER, nInstFlag::STALL, ""),
    tInstLibEntry<tMethod>("wait-cond-less", &cHardwareBCR::Inst_WaitCondition_Less, INST_CLASS_OTHER, nInstFlag::STALL, ""),
    tInstLibEntry<tMethod>("wait-cond-gtr", &cHardwareBCR::Inst_WaitCondition_Greater, INST_CLASS_OTHER, nInstFlag::STALL, ""),
        
    // State Grid instructions
    tInstLibEntry<tMethod>("sg-move", &cHardwareBCR::Inst_SGMove, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("sg-rotate-l", &cHardwareBCR::Inst_SGRotateL, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("sg-rotate-r", &cHardwareBCR::Inst_SGRotateR, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("sg-sense", &cHardwareBCR::Inst_SGSense, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_INPUT),
    
    // Movement and Navigation instructions
    tInstLibEntry<tMethod>("move", &cHardwareBCR::Inst_Move, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("get-north-offset", &cHardwareBCR::Inst_GetNorthOffset, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("get-position-offset", &cHardwareBCR::Inst_GetPositionOffset, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_INPUT),    
    tInstLibEntry<tMethod>("get-northerly", &cHardwareBCR::Inst_GetNortherly, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("get-easterly", &cHardwareBCR::Inst_GetEasterly, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("zero-easterly", &cHardwareBCR::Inst_ZeroEasterly, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("zero-northerly", &cHardwareBCR::Inst_ZeroNortherly, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("zero-position-offset", &cHardwareBCR::Inst_ZeroPosOffset, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_INPUT),

    // Rotation
    tInstLibEntry<tMethod>("rotate-home", &cHardwareBCR::Inst_RotateHome, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("rotate-to-unoccupied-cell", &cHardwareBCR::Inst_RotateUnoccupiedCell, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("rotate-x", &cHardwareBCR::Inst_RotateX, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION),
    
    // Resource and Topography Sensing
    tInstLibEntry<tMethod>("sense-resource-id", &cHardwareBCR::Inst_SenseResourceID, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT), 
    tInstLibEntry<tMethod>("sense-nest", &cHardwareBCR::Inst_SenseNest, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("sense-faced-habitat", &cHardwareBCR::Inst_SenseFacedHabitat, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("look-ahead", &cHardwareBCR::Inst_LookAhead, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("look-ahead-intercept", &cHardwareBCR::Inst_LookAheadIntercept, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("look-around", &cHardwareBCR::Inst_LookAround, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("look-around-intercept", &cHardwareBCR::Inst_LookAroundIntercept, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("look-ft", &cHardwareBCR::Inst_LookFT, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
    tInstLibEntry<tMethod>("look-around-ft", &cHardwareBCR::Inst_LookAroundFT, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
    
    tInstLibEntry<tMethod>("set-forage-target", &cHardwareBCR::Inst_SetForageTarget, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("set-ft-once", &cHardwareBCR::Inst_SetForageTargetOnce, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("set-rand-ft-once", &cHardwareBCR::Inst_SetRandForageTargetOnce, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("get-forage-target", &cHardwareBCR::Inst_GetForageTarget, INST_CLASS_ENVIRONMENT, 0, "", BEHAV_CLASS_INPUT),
    
    tInstLibEntry<tMethod>("collect-specific", &cHardwareBCR::Inst_CollectSpecific, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("get-res-stored", &cHardwareBCR::Inst_GetResStored, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
 
    // Opinion instructions.
    tInstLibEntry<tMethod>("set-opinion", &cHardwareBCR::Inst_SetOpinion, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("get-opinion", &cHardwareBCR::Inst_GetOpinion, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),

    // Grouping instructions
    tInstLibEntry<tMethod>("join-group", &cHardwareBCR::Inst_JoinGroup, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION),
    tInstLibEntry<tMethod>("get-group-id", &cHardwareBCR::Inst_GetGroupID, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
    
    // Org Interaction instructions
    tInstLibEntry<tMethod>("get-faced-org-id", &cHardwareBCR::Inst_GetFacedOrgID, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
    
    tInstLibEntry<tMethod>("teach-offspring", &cHardwareBCR::Inst_TeachOffspring, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION), 
    tInstLibEntry<tMethod>("learn-parent", &cHardwareBCR::Inst_LearnParent, INST_CLASS_ENVIRONMENT, nInstFlag::STALL, "", BEHAV_CLASS_ACTION), 

    // Control-type Instructions
    tInstLibEntry<tMethod>("scramble-registers", &cHardwareBCR::Inst_ScrambleReg, INST_CLASS_DATA, nInstFlag::STALL, "", BEHAV_CLASS_INPUT),
  };
  
  
  const int n_size = sizeof(s_n_array)/sizeof(cNOPEntry);
  
  static cString n_names[n_size];
  static int nop_mods[n_size];
  for (int i = 0; i < n_size && i < NUM_REGISTERS; i++) {
    n_names[i] = s_n_array[i].name;
    nop_mods[i] = s_n_array[i].nop_mod;
  }
  
  const int f_size = sizeof(s_f_array)/sizeof(tInstLibEntry<tMethod>);
  static tMethod functions[f_size];
  for (int i = 0; i < f_size; i++) functions[i] = s_f_array[i].GetFunction();
  
  const int def = 0;
  const int null_inst = 16;
  
  return new tInstLib<tMethod>(f_size, s_f_array, n_names, nop_mods, functions, def, null_inst);
}

cHardwareBCR::cHardwareBCR(cAvidaContext& ctx, cWorld* world, cOrganism* in_organism, cInstSet* in_inst_set)
: cHardwareBase(world, in_organism, in_inst_set), m_mem_array(1), m_sensor(world, in_organism)
{
  m_functions = s_inst_slib->GetFunctions();
  
  m_spec_die = false;
  
  m_thread_slicing_parallel = (m_world->GetConfig().THREAD_SLICING_METHOD.Get() == 1);
  m_no_cpu_cycle_time = m_world->GetConfig().NO_CPU_CYCLE_TIME.Get();
  
  m_slip_read_head = !m_world->GetConfig().SLIP_COPY_MODE.Get();
  
  const Genome& in_genome = in_organism->GetGenome();
  ConstInstructionSequencePtr in_seq_p;
  in_seq_p.DynamicCastFrom(in_genome.Representation());
  const InstructionSequence& in_seq = *in_seq_p;
  
  m_mem_array[0] = in_seq;  // Initialize memory...
  m_use_avatar = m_world->GetConfig().USE_AVATARS.Get();
  Reset(ctx);                            // Setup the rest of the hardware...
}


void cHardwareBCR::internalReset()
{
  m_cycle_count = 0;
  m_last_output = 0;
  
  m_use_avatar = m_world->GetConfig().USE_AVATARS.Get();
  m_sensor.Reset();
  
  
  // Stack
  m_global_stack.Clear(m_inst_set->GetStackSize());
  
  
  // Memory
  m_mem_array.Resize(1);
  for (int i = 1; i < MAX_MEM_SPACES; i++) m_mem_ids[i] = -1;
  m_mem_ids[0] = 0;
  
  
  // Threads
  m_threads.Resize(1);
  for (int i = 1; i < MAX_THREADS; i++) m_thread_ids[i] = -1;
  
  m_threads[0].Reset(this, 0);
  m_thread_ids[0] = 0;
  m_cur_thread = 0;
  m_waiting_threads = 0;
  m_running_threads = 1;
}


void cHardwareBCR::internalResetOnFailedDivide()
{
	internalReset();
}


void cHardwareBCR::Thread::Reset(cHardwareBCR* in_hardware, int in_id)
{
  thread_id = in_id;
  
  for (int i = 0; i < NUM_REGISTERS; i++) reg[i].Clear();
  for (int i = 0; i < NUM_HEADS; i++) heads[i].Reset(in_hardware);
  stack.Clear(in_hardware->GetInstSet().GetStackSize());
  
  cur_stack = 0;
  
  reading_label = false;
  reading_seq = false;
  running = true;
  active = true;
  read_label.Clear();
  next_label.Clear();
}


// This function processes the very next command in the genome, and is made
// to be as optimized as possible.  This is the heart of avida.

bool cHardwareBCR::SingleProcess(cAvidaContext& ctx, bool speculative)
{
  assert(!speculative || (speculative && !m_thread_slicing_parallel));
  
  // Mark this organism as running...
  m_organism->SetRunning(true);
  
  if (!speculative && m_spec_die) {
    m_organism->Die(ctx);
    m_organism->SetRunning(false);
    return false;
  }
  
  cPhenotype& phenotype = m_organism->GetPhenotype();
  
  m_cycle_count++;
  phenotype.IncCPUCyclesUsed();
  if (!m_no_cpu_cycle_time) phenotype.IncTimeUsed();


  // Wake any stalled threads
  for (int i = 0; i < m_threads.GetSize(); i++) {
    if (!m_threads[i].active && m_threads[i].wait_reg == -1) m_threads[i].active = true;
  }

  bool speculative_stall = false;
  bool behav_class_used[3] = { false, false, false };
  
  // Execute specified number of micro ops per cpu cycle on each thread in a round robin fashion
  const int uop_ratio = m_inst_set->GetUOpsPerCycle();
  for (int i = 0; i < uop_ratio; i++) {
    for (m_cur_thread = 0; m_cur_thread < m_threads.GetSize(); m_cur_thread++) {
      // Setup the hardware for the next instruction to be executed.
      
      // If the currently selected thread is inactive, proceed to the next thread
      if (!m_threads[m_cur_thread].running || !m_threads[m_cur_thread].active) continue;
      
      m_advance_ip = true;
      cHeadCPU& ip = m_threads[m_cur_thread].heads[hIP];
      ip.Adjust();
    
      // Print the status of this CPU at each step...
      if (m_tracer != NULL) m_tracer->TraceHardware(ctx, *this);
      
      // Find the instruction to be executed
      const Instruction& cur_inst = ip.GetInst();
      
      if (speculative && (m_spec_die || m_inst_set->ShouldStall(cur_inst))) {
        // Speculative instruction stall, flag it and halt the thread
        speculative_stall = true;
        m_threads[m_cur_thread].active = false;
        m_threads[m_cur_thread].wait_reg = -1;
        continue;
      }
      
      // Print the short form status of this CPU at each step...
      if (m_minitracer != NULL) m_minitracer->TraceHardware(ctx, *this, false, true);
      
      bool exec = true;
      int exec_success = 0;

      BehavClass behav_class = m_inst_set->GetInstLib()->Get(m_inst_set->GetLibFunctionIndex(ip.GetInst())).GetBehavClass();
      
      // Check if this instruction class has been used and should cause the thread to stall?
      if (behav_class < BEHAV_CLASS_NONE && behav_class_used[behav_class]) {
        m_threads[m_cur_thread].active = false;
        m_threads[m_cur_thread].wait_reg = -1;
        continue;
      }

      // Test if costs have been paid and it is okay to execute this now...
      
      // record any failure due to costs being paid
      // before we try to execute the instruction, is this org currently paying precosts for it
      bool on_pause = IsPayingActiveCost(ctx, m_cur_thread);
      if (m_has_any_costs) exec = SingleProcess_PayPreCosts(ctx, cur_inst, m_cur_thread);
      if (!exec) exec_success = -1;
      
      // Now execute the instruction...
      bool rand_fail = false;
      if (exec == true) {
        // NOTE: This call based on the cur_inst must occur prior to instruction
        //       execution, because this instruction reference may be invalid after
        //       certain classes of instructions (namely divide instructions) @DMB
        const int addl_time_cost = m_inst_set->GetAddlTimeCost(cur_inst);
        
        // Prob of exec (moved from SingleProcess_PayCosts so that we advance IP after a fail)
        if ( m_inst_set->GetProbFail(cur_inst) > 0.0 ) {
          exec = !( ctx.GetRandom().P(m_inst_set->GetProbFail(cur_inst)) );
          rand_fail = !exec;
        }
        
        if (exec == true) {
          if (SingleProcess_ExecuteInst(ctx, cur_inst)) {
            SingleProcess_PayPostResCosts(ctx, cur_inst);
            SingleProcess_SetPostCPUCosts(ctx, cur_inst, m_cur_thread);
            // record execution success
            exec_success = 1;
          }
        }
        
        // Check if the instruction just executed caused premature death, break out of execution if so
        if (phenotype.GetToDelete()) {
          if (m_minitracer != NULL) m_minitracer->TraceHardware(ctx, *this, false, true, exec_success);
          break;
        }
        
        // Some instruction (such as jump) may turn m_advance_ip off.  Usually
        // we now want to move to the next instruction in the memory.
        if (m_advance_ip == true) ip.Advance();
        
        // Pay the additional death_cost of the instruction now
        phenotype.IncTimeUsed(addl_time_cost);


        // mark behavior class as used, when appropriate
        if (behav_class < BEHAV_CLASS_NONE) behav_class_used[behav_class] = true;
      }
      
      // if using mini traces, report success or failure of execution
      if (m_minitracer != NULL) m_minitracer->TraceHardware(ctx, *this, false, true, exec_success);
      
      bool do_record = false;
      // record exec failed if the org just now started paying precosts
      if (exec_success == -1 && !on_pause) do_record = true;
      // if exec succeeded but was on pause before this execution, we already recorded it
      // otherwise we record what the org did
      else if (exec_success == 1 && !on_pause) do_record = true;
      // if random failure, we record 'what the org was trying to do'
      else if (rand_fail) do_record = true;
      // if exec failed because of something inside the instruction itself, record the attempt
      else if (exec_success == 0) do_record = true;
      if (do_record) {
        // this will differ from time used
        phenotype.IncNumExecs();
        if (m_microtrace || m_topnavtrace) {
          RecordMicroTrace(cur_inst);
          if (m_topnavtrace) RecordNavTrace(m_use_avatar);
        }
      }
      
      if (phenotype.GetToDelete()) break;
    }
    
    if (phenotype.GetToDelete()) break;
  }

  // Kill creatures who have reached their max num of instructions executed
  const int max_executed = m_organism->GetMaxExecuted();
  if ((max_executed > 0 && phenotype.GetTimeUsed() >= max_executed) || phenotype.GetToDie() == true) {
    if (speculative) m_spec_die = true;
    else m_organism->Die(ctx);
  }
  if (!speculative && phenotype.GetToDelete()) m_spec_die = true;
  
  m_organism->SetRunning(false);
  CheckImplicitRepro(ctx);
  
  return !m_spec_die && !speculative_stall;
}

// This method will handle the actuall execution of an instruction
// within single process, once that function has been finalized.
bool cHardwareBCR::SingleProcess_ExecuteInst(cAvidaContext& ctx, const Instruction& cur_inst) 
{
  // Copy Instruction locally to handle stochastic effects
  Instruction actual_inst = cur_inst;
  
#ifdef EXECUTION_ERRORS
  // If there is an execution error, execute a random instruction.
  if (m_organism->TestExeErr()) actual_inst = m_inst_set->GetRandomInst(ctx);
#endif /* EXECUTION_ERRORS */
  
  // Get a pointer to the corrisponding method...
  int inst_idx = m_inst_set->GetLibFunctionIndex(actual_inst);
  
  // Mark the instruction as executed
  getIP().SetFlagExecuted();
	
  // instruction execution count incremeneted
  m_organism->GetPhenotype().IncCurInstCount(actual_inst.GetOp());
  
  // And execute it.
  const bool exec_success = (this->*(m_functions[inst_idx]))(ctx);
  
  // decremenet if the instruction was not executed successfully
  if (exec_success == false) {
    m_organism->GetPhenotype().DecCurInstCount(actual_inst.GetOp());
  }  
  return exec_success;
}


void cHardwareBCR::ProcessBonusInst(cAvidaContext& ctx, const Instruction& inst)
{
  // Mark this organism as running...
  bool prev_run_state = m_organism->IsRunning();
  m_organism->SetRunning(true);
  
  if (m_tracer != NULL) m_tracer->TraceHardware(ctx, *this, true);
  
  SingleProcess_ExecuteInst(ctx, inst);
  
  m_organism->SetRunning(prev_run_state);
}


void cHardwareBCR::PrintStatus(ostream& fp)
{
  fp << "CPU CYCLE:" << m_organism->GetPhenotype().GetCPUCyclesUsed() << " ";
  fp << "THREAD:" << m_cur_thread << "  ";
  fp << "IP:" << getIP().GetPosition() << "    ";
  
  
  for (int i = 0; i < NUM_REGISTERS; i++) {
    DataValue& reg = m_threads[m_cur_thread].reg[i];
    fp << static_cast<char>('A' + i) << "X:" << GetRegister(i) << " ";
    fp << setbase(16) << "[0x" << reg.value <<  "] " << setbase(10);
    fp << "(" << reg.from_env << " " << reg.env_component << " " << reg.originated << " " << reg.oldest_component << ")  ";
  }
  
  // Add some extra information if additional time costs are used for instructions,
  // leave this out if there are no differences to keep it cleaner
  if ( m_organism->GetPhenotype().GetTimeUsed() != m_organism->GetPhenotype().GetCPUCyclesUsed() )
  {
    fp << "  EnergyUsed:" << m_organism->GetPhenotype().GetTimeUsed();
  }
  fp << endl;
  
  fp << "  R-Head:" << getHead(hREAD).GetPosition() << " "
  << "W-Head:" << getHead(hWRITE).GetPosition()  << " "
  << "F-Head:" << getHead(hFLOW).GetPosition()   << "  "
  << "RL:" << GetReadLabel().AsString() << "   "
  << "Ex:" << m_last_output
  << endl;
  
  int number_of_stacks = GetNumStacks();
  for (int stack_id = 0; stack_id < number_of_stacks; stack_id++) {
    fp << ((m_threads[m_cur_thread].cur_stack == stack_id) ? '*' : ' ') << " Stack " << stack_id << ":" << setbase(16) << setfill('0');
    for (int i = 0; i < m_inst_set->GetStackSize(); i++) fp << " Ox" << setw(8) << GetStack(i, stack_id, 0);
    fp << setfill(' ') << setbase(10) << endl;
  }
  
  for (int i = 0; i < m_mem_array.GetSize(); i++) {
    const cCPUMemory& mem = m_mem_array[i];
    fp << "  Mem " << i << " (" << mem.GetSize() << "): " << mem.AsString() << endl;
  }
  
  fp.flush();
}

void cHardwareBCR::SetupMiniTraceFileHeader(const cString& filename, const int gen_id, const cString& genotype)
{
  const Genome& in_genome = m_organism->GetGenome();
  ConstInstructionSequencePtr in_seq_p;
  in_seq_p.DynamicCastFrom(in_genome.Representation());
  const InstructionSequence& in_seq = *in_seq_p;

  cDataFile& df = m_world->GetDataFile(filename);
  df.WriteTimeStamp();
  cString org_dat("");
  df.WriteComment(org_dat.Set("Update Born: %d", m_world->GetStats().GetUpdate()));
  df.WriteComment(org_dat.Set("Org ID: %d", m_organism->GetID()));
  df.WriteComment(org_dat.Set("Genotype ID: %d", gen_id));
  df.WriteComment(org_dat.Set("Genotype: %s", (const char*) genotype));
  df.WriteComment(org_dat.Set("Genome Length: %d", in_seq.GetSize()));
  df.WriteComment(" ");
  df.WriteComment("Exec Stats Columns:");
  df.WriteComment("CPU Cycle");
  df.WriteComment("Current Update");
  df.WriteComment("Register Contents (CPU Cycle Origin of Contents)");
  df.WriteComment("Current Thread");
  df.WriteComment("IP Position");
  df.WriteComment("RH Position");
  df.WriteComment("WH Position");
  df.WriteComment("FH Position");
  df.WriteComment("CPU Cycle of Last Output");
  df.WriteComment("Current Merit");
  df.WriteComment("Current Bonus");
  df.WriteComment("Forager Type");
  df.WriteComment("Group ID (opinion)");
  df.WriteComment("Current Cell");
  df.WriteComment("Avatar Cell");
  df.WriteComment("Faced Direction");
  df.WriteComment("Faced Cell Occupied?");
  df.WriteComment("Faced Cell Has Hill?");
  df.WriteComment("Faced Cell Has Wall?");
  df.WriteComment("Queued Instruction");
  df.WriteComment("Trailing NOPs");
  df.WriteComment("Did Queued Instruction Execute (-1=no, paying cpu costs; 0=failed; 1=yes)");
  df.Endl();
}

void cHardwareBCR::PrintMiniTraceStatus(cAvidaContext& ctx, ostream& fp, const cString& next_name)
{
  // basic status info
  fp << m_cycle_count << " ";
  fp << m_world->GetStats().GetUpdate() << " ";
  for (int i = 0; i < NUM_REGISTERS; i++) {
    DataValue& reg = m_threads[m_cur_thread].reg[i];
    fp << GetRegister(i) << " ";
    fp << "(" << reg.originated << ") ";
  }    
  // genome loc info
  fp << m_cur_thread << " ";
  fp << getIP().GetPosition() << " ";  
  fp << getHead(hREAD).GetPosition() << " ";
  fp << getHead(hWRITE).GetPosition()  << " ";
  fp << getHead(hFLOW).GetPosition()   << " ";
  // last output
  fp << m_last_output << " ";
  // phenotype/org status info
  fp << m_organism->GetPhenotype().GetMerit().GetDouble() << " ";
  fp << m_organism->GetPhenotype().GetCurBonus() << " ";
  fp << m_organism->GetForageTarget() << " ";
  if (m_organism->HasOpinion()) fp << m_organism->GetOpinion().first << " ";
  else fp << -99 << " ";
  // environment info / things that affect movement
  fp << m_organism->GetOrgInterface().GetCellID() << " ";
  if (m_use_avatar) fp << m_organism->GetOrgInterface().GetAVCellID() << " ";
  if (!m_use_avatar) fp << m_organism->GetOrgInterface().GetFacedDir() << " ";
  else fp << m_organism->GetOrgInterface().GetAVFacing() << " ";
  if (!m_use_avatar) fp << m_organism->IsNeighborCellOccupied() << " ";  
  else fp << m_organism->GetOrgInterface().FacedHasAV() << " ";
  const cResourceLib& resource_lib = m_world->GetEnvironment().GetResourceLib();
  Apto::Array<double> cell_resource_levels;
  if (!m_use_avatar) cell_resource_levels = m_organism->GetOrgInterface().GetFacedCellResources(ctx);
  else cell_resource_levels = m_organism->GetOrgInterface().GetAVFacedResources(ctx);
  int wall = 0;
  int hill = 0;
  for (int i = 0; i < cell_resource_levels.GetSize(); i++) {
    if (resource_lib.GetResource(i)->GetHabitat() == 2 && cell_resource_levels[i] > 0) wall = 1;
    if (resource_lib.GetResource(i)->GetHabitat() == 1 && cell_resource_levels[i] > 0) hill = 1;
    if (hill == 1 && wall == 1) break;
  }
  fp << hill << " ";
  fp << wall << " ";
  // instruction about to be executed
  fp << next_name << " ";
  // any trailing nops (up to NUM_REGISTERS)
  cCPUMemory& memory = m_mem_array[0];
  int pos = getIP().GetPosition();
  Apto::Array<int, Apto::Smart> seq;
  seq.Resize(0);
  for (int i = 0; i < NUM_REGISTERS; i++) {
    pos += 1;
    if (pos >= memory.GetSize()) pos = 0;
    if (m_inst_set->IsNop(memory[pos])) seq.Push(m_inst_set->GetNopMod(memory[pos])); 
    else break;
  }
  cString mod_string;
  for (int j = 0; j < seq.GetSize(); j++) {
    mod_string += (char) seq[j] + 'A';  
  }  
  if (mod_string.GetSize() != 0) fp << mod_string << " ";
  else fp << "NoMods" << " ";
}

void cHardwareBCR::PrintMiniTraceSuccess(ostream& fp, const int exec_sucess)
{
  fp << exec_sucess;
  fp << endl;
  fp.flush();
}

void cHardwareBCR::FindLabelStart(cHeadCPU& head, bool mark_executed)
{
  cHeadCPU& ip = getIP();
  const cCodeLabel& search_label = GetLabel();
  
  // Make sure the label is of size > 0.
  if (search_label.GetSize() == 0) {
    head.Set(ip);
    return;
  }
  
  cCPUMemory& memory = head.GetMemory();
  int pos = 0;
  
  while (pos < memory.GetSize()) {
    if (m_inst_set->IsLabel(memory[pos])) { // starting label found
      pos++;
      
      // Check for direct matched label pattern, can be substring of 'label'ed target
      // - must match all NOPs in search_label
      // - extra NOPs in 'label'ed target are ignored
      int size_matched = 0;
      while (size_matched < search_label.GetSize() && pos < memory.GetSize()) {
        if (!m_inst_set->IsNop(memory[pos]) || search_label[size_matched] != m_inst_set->GetNopMod(memory[pos])) break;
        size_matched++;
        pos++;
      }
      
      // Check that the label matches and has examined the full sequence of nops following the 'label' instruction
      if (size_matched == search_label.GetSize()) {
        // Return Head pointed at last NOP of label sequence
        if (mark_executed) {
          size_matched++; // Increment size matched so that it includes the label instruction
          const int start = pos - size_matched;
          const int max = m_world->GetConfig().MAX_LABEL_EXE_SIZE.Get() + 1; // Max label + 1 for the label instruction itself
          for (int i = 0; i < size_matched && i < max; i++) memory.SetFlagExecuted(start + i);
        }
        head.Set(pos - 1, head.GetMemSpace());
        return;
      }
      
      continue; 
    }
    pos++;
  }
  
  // Return start point if not found
  head.Set(ip);
}

void cHardwareBCR::FindNopSequenceStart(cHeadCPU& head, bool mark_executed)
{
  cHeadCPU& ip = getIP();
  const cCodeLabel& search_label = GetLabel();
  
  // Make sure the label is of size > 0.
  if (search_label.GetSize() == 0) {
    head.Set(ip);
    return;
  }
  
  cCPUMemory& memory = head.GetMemory();
  int pos = 0;
  
  while (pos < memory.GetSize()) {
    if (m_inst_set->IsNop(memory[pos])) { // start of sequence found
      
      // Check for direct matched label pattern, can be substring of 'label'ed target
      // - must match all NOPs in search_label
      // - extra NOPs in 'label'ed target are ignored
      int size_matched = 0;
      while (size_matched < search_label.GetSize() && pos < memory.GetSize()) {
        if (!m_inst_set->IsNop(memory[pos]) || search_label[size_matched] != m_inst_set->GetNopMod(memory[pos])) break;
        size_matched++;
        pos++;
      }
      
      // Check that the label matches and has examined the full sequence of nops following the 'label' instruction
      if (size_matched == search_label.GetSize()) {
        // Return Head pointed at last NOP of label sequence
        if (mark_executed) {
          const int start = pos - size_matched;
          const int max = m_world->GetConfig().MAX_LABEL_EXE_SIZE.Get();
          for (int i = 0; i < size_matched && i < max; i++) memory.SetFlagExecuted(start + i);
        }
        head.Set(pos - 1, head.GetMemSpace());
        return;
      }
    }
    pos++;
  }
  
  // Return start point if not found
  head.Set(ip);
}


void cHardwareBCR::FindLabelForward(cHeadCPU& head, bool mark_executed)
{
  cHeadCPU& ip = getIP();
  const cCodeLabel& search_label = GetLabel();
  
  // Make sure the label is of size > 0.
  if (search_label.GetSize() == 0) {
    head.Set(ip);
    return;
  }
  
  cHeadCPU pos(head);
  pos++;
  
  while (pos.GetPosition() != head.GetPosition()) {
    if (m_inst_set->IsLabel(pos.GetInst())) { // starting label found
      const int label_start = pos.GetPosition();
      pos++;
      
      // Check for direct matched label pattern, can be substring of 'label'ed target
      // - must match all NOPs in search_label
      // - extra NOPs in 'label'ed target are ignored
      int size_matched = 0;
      while (size_matched < search_label.GetSize() && pos.GetPosition() != head.GetPosition()) {
        if (!m_inst_set->IsNop(pos.GetInst()) || search_label[size_matched] != m_inst_set->GetNopMod(pos.GetInst())) break;
        size_matched++;
        pos++;
      }
      
      // Check that the label matches and has examined the full sequence of nops following the 'label' instruction
      if (size_matched == search_label.GetSize()) {
        pos--;
        const int found_pos = pos.GetPosition();
        
        if (mark_executed) {
          pos.Set(label_start);
          const int max = m_world->GetConfig().MAX_LABEL_EXE_SIZE.Get() + 1; // Max label + 1 for the label instruction itself
          for (int i = 0; i < size_matched && i < max; i++, pos++) pos.SetFlagExecuted();
        }
        
        // Return Head pointed at last NOP of label sequence
        head.Set(found_pos, head.GetMemSpace());
        return;
      }
      
      continue; 
    }
    pos++;
  }
  
  // Return start point if not found
  head.Set(ip);
}

void cHardwareBCR::FindLabelBackward(cHeadCPU& head, bool mark_executed)
{
  cHeadCPU& ip = getIP();
  const cCodeLabel& search_label = GetLabel();
  
  // Make sure the label is of size > 0.
  if (search_label.GetSize() == 0) {
    head.Set(ip);
    return;
  }
  
  cHeadCPU lpos(head);
  cHeadCPU pos(head);
  lpos--;
  
  while (pos.GetPosition() != head.GetPosition()) {
    if (m_inst_set->IsLabel(lpos.GetInst())) { // starting label found
      pos.Set(lpos.GetPosition());
      pos++;
      
      // Check for direct matched label pattern, can be substring of 'label'ed target
      // - must match all NOPs in search_label
      // - extra NOPs in 'label'ed target are ignored
      int size_matched = 0;
      while (size_matched < search_label.GetSize() && pos.GetPosition() != head.GetPosition()) {
        if (!m_inst_set->IsNop(pos.GetInst()) || search_label[size_matched] != m_inst_set->GetNopMod(pos.GetInst())) break;
        size_matched++;
        pos++;
      }
      
      // Check that the label matches and has examined the full sequence of nops following the 'label' instruction
      if (size_matched == search_label.GetSize()) {
        pos--;
        const int found_pos = pos.GetPosition();
        
        if (mark_executed) {
          const int max = m_world->GetConfig().MAX_LABEL_EXE_SIZE.Get() + 1; // Max label + 1 for the label instruction itself
          for (int i = 0; i < size_matched && i < max; i++, lpos++) lpos.SetFlagExecuted();
        }
        
        // Return Head pointed at last NOP of label sequence
        head.Set(found_pos, head.GetMemSpace());
        return;
      }
    }
    lpos--;
  }
  
  // Return start point if not found
  head.Set(ip);
}




void cHardwareBCR::FindNopSequenceForward(cHeadCPU& head, bool mark_executed)
{
  cHeadCPU& ip = getIP();
  const cCodeLabel& search_label = GetLabel();
  
  // Make sure the label is of size > 0.
  if (search_label.GetSize() == 0) {
    head.Set(ip);
    return;
  }
  
  cHeadCPU pos(head);
  pos++;
  
  while (pos.GetPosition() != head.GetPosition()) {
    if (m_inst_set->IsNop(pos.GetInst())) { // starting label found
      const int label_start = pos.GetPosition();
      
      // Check for direct matched nop sequence, can be substring of target
      // - must match all NOPs in search_label
      // - extra NOPs in target are ignored
      int size_matched = 0;
      while (size_matched < search_label.GetSize() && pos.GetPosition() != head.GetPosition()) {
        if (!m_inst_set->IsNop(pos.GetInst()) || search_label[size_matched] != m_inst_set->GetNopMod(pos.GetInst())) break;
        size_matched++;
        pos++;
      }
      
      // Check that the label matches and has examined the full sequence of nops
      if (size_matched == search_label.GetSize()) {
        pos--;
        const int found_pos = pos.GetPosition();
        
        if (mark_executed) {
          pos.Set(label_start);
          const int max = m_world->GetConfig().MAX_LABEL_EXE_SIZE.Get();
          for (int i = 0; i < size_matched && i < max; i++, pos++) pos.SetFlagExecuted();
        }
        
        // Return Head pointed at last NOP of label sequence
        head.Set(found_pos, head.GetMemSpace());
        return;
      }
    }
    
    if (pos.GetPosition() == head.GetPosition()) break;
    pos++;
  }
  
  // Return start point if not found
  head.Set(ip);
}


void cHardwareBCR::FindNopSequenceBackward(cHeadCPU& head, bool mark_executed)
{
  cHeadCPU& ip = getIP();
  const cCodeLabel& search_label = GetLabel();
  
  // Make sure the label is of size > 0.
  if (search_label.GetSize() == 0) {
    head.Set(ip);
    return;
  }
  
  cHeadCPU lpos(head);
  cHeadCPU pos(head);
  lpos--;
  
  while (pos.GetPosition() != head.GetPosition()) {
    if (m_inst_set->IsNop(pos.GetInst())) { // starting label found
      pos.Set(lpos.GetPosition());
      
      // Check for direct matched nop sequence, can be substring of target
      // - must match all NOPs in search_label
      // - extra NOPs in target are ignored
      int size_matched = 0;
      while (size_matched < search_label.GetSize() && pos.GetPosition() != head.GetPosition()) {
        if (!m_inst_set->IsNop(pos.GetInst()) || search_label[size_matched] != m_inst_set->GetNopMod(pos.GetInst())) break;
        size_matched++;
        pos++;
      }
      
      // Check that the label matches and has examined the full sequence of nops
      if (size_matched == search_label.GetSize()) {
        pos--;
        const int found_pos = pos.GetPosition();
        
        if (mark_executed) {
          const int max = m_world->GetConfig().MAX_LABEL_EXE_SIZE.Get();
          for (int i = 0; i < size_matched && i < max; i++, lpos++) lpos.SetFlagExecuted();
        }
        
        // Return Head pointed at last NOP of label sequence
        head.Set(found_pos, head.GetMemSpace());
        return;
      }
      
      continue; 
    }
    lpos--;
  }
  
  // Return start point if not found
  head.Set(ip);
}

void cHardwareBCR::ReadInst(Instruction in_inst)
{
  bool is_nop = m_inst_set->IsNop(in_inst);
  
  if (m_inst_set->IsLabel(in_inst)) {
    GetReadLabel().Clear();
    m_threads[m_cur_thread].reading_label = true;
  } else if (m_threads[m_cur_thread].reading_label && is_nop) {
    GetReadLabel().AddNop(in_inst.GetOp());
  } else {
    GetReadLabel().Clear();
    m_threads[m_cur_thread].reading_label = false;
  }
  
  if (!m_threads[m_cur_thread].reading_seq && is_nop) {
    GetReadSequence().AddNop(in_inst.GetOp());
    m_threads[m_cur_thread].reading_seq = true;
  } else if (m_threads[m_cur_thread].reading_seq && is_nop) {
    GetReadSequence().AddNop(in_inst.GetOp());
  } else {
    GetReadSequence().Clear();
    m_threads[m_cur_thread].reading_seq = false;
  }
}

// This function looks at the current position in the info of a creature,
// and sets the next_label to be the sequence of nops which follows.  The
// instruction pointer is left on the last line of the label found.

void cHardwareBCR::ReadLabel(int max_size)
{
  int count = 0;
  cHeadCPU * inst_ptr = &( getIP() );
  
  GetLabel().Clear();
  
  while (m_inst_set->IsNop(inst_ptr->GetNextInst()) &&
         (count < max_size)) {
    count++;
    inst_ptr->Advance();
    GetLabel().AddNop(m_inst_set->GetNopMod(inst_ptr->GetInst()));
    
    // If this is the first line of the template, mark it executed.
    if (GetLabel().GetSize() <=	m_world->GetConfig().MAX_LABEL_EXE_SIZE.Get()) {
      inst_ptr->SetFlagExecuted();
    }
  }
}

bool cHardwareBCR::ThreadCreate(int thread_label, const cHeadCPU& start_pos)
{
  
  int thread_id = m_threads.GetSize();
  
  // Check for thread cap, base thread label (i.e. no label)
  if (thread_id == m_world->GetConfig().MAX_CPU_THREADS.Get()) return false;
  
  // Check for existing thread
  if (thread_label >= 0 && m_thread_ids[thread_label] >= 0) {
    thread_id = m_thread_ids[thread_label];
    if (m_threads[thread_id].running) {
      return -1;  // Thread exists, and is running... call fails
    } else {
      m_threads[thread_id].Reset(this, thread_id);
      m_threads[thread_id].heads[hIP] = start_pos;
      m_running_threads++;
      return thread_id;
    }
  }
  
  // Add new thread entry
  m_threads.Resize(thread_id + 1);

  if (thread_label >= 0) {
    m_thread_ids[thread_label] = thread_id;
  } else {
    for (int i = 0; i < MAX_THREADS; i++) {
      if (m_thread_ids[i] < 0) {
        m_thread_ids[i] = thread_id;
        break;
      }
    }
  }
  
  // Setup this thread into the current selected memory space (Flow Head)
  m_threads[thread_id].Reset(this, thread_id);
  m_threads[thread_id].heads[hIP] = start_pos;
  m_running_threads++;
  return thread_id;
}


////////////////////////////
//  Instruction Helpers...
////////////////////////////

inline int cHardwareBCR::FindModifiedRegister(int default_register)
{
  assert(default_register < NUM_REGISTERS);  // Reg ID too high.
  
  if (m_inst_set->IsNop(getIP().GetNextInst())) {
    getIP().Advance();
    default_register = m_inst_set->GetNopMod(getIP().GetInst());
    getIP().SetFlagExecuted();
  }
  return default_register;
}

inline int cHardwareBCR::FindModifiedNextRegister(int default_register)
{
  assert(default_register < NUM_REGISTERS);  // Reg ID too high.
  
  if (m_inst_set->IsNop(getIP().GetNextInst())) {
    getIP().Advance();
    default_register = m_inst_set->GetNopMod(getIP().GetInst());
    getIP().SetFlagExecuted();
  } else {
    default_register = (default_register + 1) % NUM_REGISTERS;
  }
  return default_register;
}

inline int cHardwareBCR::FindModifiedPreviousRegister(int default_register)
{
  assert(default_register < NUM_REGISTERS);  // Reg ID too high.
  
  if (m_inst_set->IsNop(getIP().GetNextInst())) {
    getIP().Advance();
    default_register = m_inst_set->GetNopMod(getIP().GetInst());
    getIP().SetFlagExecuted();
  } else {
    default_register = (default_register + NUM_REGISTERS - 1) % NUM_REGISTERS;
  }
  return default_register;
}


inline int cHardwareBCR::FindModifiedHead(int default_head)
{
  assert(default_head < NUM_HEADS); // Head ID too high.
  
  if (m_inst_set->IsNop(getIP().GetNextInst())) {
    getIP().Advance();
    default_head = m_inst_set->GetNopMod(getIP().GetInst());
    getIP().SetFlagExecuted();
  }
  return default_head;
}


inline int cHardwareBCR::FindNextRegister(int base_reg)
{
  return (base_reg + 1) % NUM_REGISTERS;
}



int cHardwareBCR::calcCopiedSize(const int parent_size, const int child_size)
{
  int copied_size = 0;
  const cCPUMemory& memory = m_mem_array[m_cur_offspring];
  for (int i = 0; i < memory.GetSize(); i++) {
    if (memory.FlagCopied(i)) copied_size++;
	}
  return copied_size;
}

bool cHardwareBCR::Divide_Main(cAvidaContext& ctx, int mem_space_used, int write_head_pos, double mut_multiplier)
{  
  // Make sure the memory space we're using exists
  if (m_mem_array.GetSize() <= mem_space_used) return false;
  
  // Make sure this divide will produce a viable offspring.
  m_cur_offspring = mem_space_used; // save current child memory space for use by dependent functions (e.g. calcCopiedSize())
  if (!Divide_CheckViable(ctx, m_mem_array[0].GetSize(), write_head_pos)) return false;
  
  // Since the divide will now succeed, set up the information to be sent to the new organism
  m_mem_array[mem_space_used].Resize(write_head_pos);
  
  InstructionSequencePtr offspring_seq(new InstructionSequence(m_mem_array[mem_space_used]));
  HashPropertyMap props;
  cHardwareManager::SetupPropertyMap(props, (const char*)m_inst_set->GetInstSetName());
  Genome offspring(GetType(), props, offspring_seq);
  
  m_organism->OffspringGenome() = offspring;
	
  // Handle Divide Mutations...
  Divide_DoMutations(ctx, mut_multiplier);
	
  // Many tests will require us to run the offspring through a test CPU;
  // this is, for example, to see if mutations need to be reverted or if
  // lineages need to be updated.
  Divide_TestFitnessMeasures(ctx);
	
  // reset first time instruction costs
  for (int i = 0; i < m_inst_ft_cost.GetSize(); i++) {
    m_inst_ft_cost[i] = m_inst_set->GetFTCost(Instruction(i));
  }
	
  //bool parent_alive = m_organism->ActivateDivide(ctx);
  bool parent_alive = m_organism->ActivateDivide(ctx);
  //reset the memory of the memory space that has been divided off
  m_mem_array[mem_space_used] = InstructionSequence("a");
	
  // Division Methods:
  // 0 - DIVIDE_METHOD_OFFSPRING - Create a child, leave parent state untouched.
  // 1 - DIVIDE_METHOD_SPLIT - Create a child, completely reset state of parent.
  // 2 - DIVIDE_METHOD_BIRTH - Create a child, reset state of parent's current thread.
  
  if (parent_alive) { // If the parent is no longer alive, all of this is moot
    switch (m_world->GetConfig().DIVIDE_METHOD.Get()) {
      case DIVIDE_METHOD_SPLIT:
        Reset(ctx);  // This will wipe out all parasites on a divide.
        break;
        
      case DIVIDE_METHOD_BIRTH:
        // Reset only the calling thread's state
        for(int x = 0; x < NUM_HEADS; x++) GetHead(x).Reset(this, 0);
        for(int x = 0; x < NUM_REGISTERS; x++) setInternalValue(x, 0, false);
        if(m_world->GetConfig().INHERIT_MERIT.Get() == 0) {
          m_organism->GetPhenotype().ResetMerit();
        }
        break;
        
      case DIVIDE_METHOD_OFFSPRING:
      default:
        break;
		}
		m_advance_ip = false;
	}
	
  return true;
}

void cHardwareBCR::checkWaitingThreads(int cur_thread, int reg_num)
{
  for (int i = 0; i < m_threads.GetSize(); i++) {
    if (i != cur_thread && !m_threads[i].active && int(m_threads[i].wait_reg) == reg_num) {
      int wait_value = m_threads[i].wait_value;
      int check_value = m_threads[cur_thread].reg[reg_num].value;
      if ((m_threads[i].wait_greater && check_value > wait_value) ||
          (m_threads[i].wait_equal && check_value == wait_value) ||
          (m_threads[i].wait_less && check_value < wait_value)) {
        
        // Wake up the thread with matched condition
        m_threads[i].active = true;
        m_waiting_threads--;
        
        // Set destination register to be the check value
        DataValue& dest = m_threads[i].reg[m_threads[i].wait_dst];
        dest.value = check_value;
        dest.from_env = false;
        dest.originated = m_cycle_count;
        dest.oldest_component = m_threads[cur_thread].reg[reg_num].oldest_component;
        dest.env_component = m_threads[cur_thread].reg[reg_num].env_component;
        
        // Cascade check
        if (m_waiting_threads) checkWaitingThreads(i, m_threads[i].wait_dst);
      }
    }
  }
}

//////////////////////////
// And the instructions...
//////////////////////////


// Multi-threading.
bool cHardwareBCR::Inst_ThreadCreate(cAvidaContext&)
{
  int thread_label = FindModifiedRegister(-1);
  int head_used = FindModifiedHead(hFLOW);
  
  bool success = ThreadCreate(thread_label, m_threads[m_cur_thread].heads[head_used]);
  if (!success) m_organism->Fault(FAULT_LOC_THREAD_FORK, FAULT_TYPE_FORK_TH);
  
  if (success) m_organism->GetPhenotype().SetIsMultiThread();
  
  return success;
}

bool cHardwareBCR::Inst_ThreadCancel(cAvidaContext& ctx)
{
  if (m_running_threads > 1) {
    m_threads[m_cur_thread].running = false;
    m_running_threads--;
  }
  return true;
}

bool cHardwareBCR::Inst_ThreadID(cAvidaContext&)
{
  const int reg_used = FindModifiedRegister(rBX);
  setInternalValue(reg_used, GetCurThreadID(), false);
  return true;
}

bool cHardwareBCR::Inst_Yield(cAvidaContext&)
{
  m_threads[m_cur_thread].active = false;
  m_threads[m_cur_thread].wait_reg = -1;
  return true;
}

bool cHardwareBCR::Inst_Label(cAvidaContext&)
{
  ReadLabel();
  return true;
}

bool cHardwareBCR::Inst_IfNEqu(cAvidaContext&) // Execute next if bx != ?cx?
{
  const int op1 = FindModifiedRegister(rBX);
  const int op2 = FindModifiedNextRegister(op1);
  if (GetRegister(op1) == GetRegister(op2))  getIP().Advance();
  return true;
}

bool cHardwareBCR::Inst_IfLess(cAvidaContext&) // Execute next if ?bx? < ?cx?
{
  const int op1 = FindModifiedRegister(rBX);
  const int op2 = FindModifiedNextRegister(op1);
  if (GetRegister(op1) >=  GetRegister(op2))  getIP().Advance();
  return true;
}

bool cHardwareBCR::Inst_IfNotZero(cAvidaContext&)  // Execute next if ?bx? != 0
{
  const int op1 = FindModifiedRegister(rBX);
  if (GetRegister(op1) == 0)  getIP().Advance();
  return true;
}
bool cHardwareBCR::Inst_IfEqualZero(cAvidaContext&)  // Execute next if ?bx? == 0
{
  const int op1 = FindModifiedRegister(rBX);
  if (GetRegister(op1) != 0)  getIP().Advance();
  return true;
}
bool cHardwareBCR::Inst_IfGreaterThanZero(cAvidaContext&)  // Execute next if ?bx? > 0
{
  const int op1 = FindModifiedRegister(rBX);
  if (GetRegister(op1) <= 0)  getIP().Advance();
  return true;
}

bool cHardwareBCR::Inst_IfLessThanZero(cAvidaContext&)  // Execute next if ?bx? < 0
{
  const int op1 = FindModifiedRegister(rBX);
  if (GetRegister(op1) >= 0)  getIP().Advance();
  return true;
}


bool cHardwareBCR::Inst_IfGtrX(cAvidaContext&)       // Execute next if BX > X; X value set according to NOP label
{
  // Compares value in BX to a specific value.  The value to compare to is determined by the nop label as follows:
  //    no nop label (default): valueToCompare = 1;
  //    nop-A: toggles valueToCompare sign-bit 
  //    nop-B: valueToCompare left-shift by 1-bit
  //    nop-C: valueToCompare left-shift by 2-bits
  //    nop-D: valueToCompare left-shift by 3-bits, etc.
  
  int valueToCompare = 1;
  
  ReadLabel();
  const cCodeLabel& shift_label = GetLabel();
  for (int i = 0; i < shift_label.GetSize(); i++) {
    if (shift_label[i] == rAX) {
      valueToCompare *= -1;
    } else {
      valueToCompare <<= shift_label[i];
    }
  }
  
  if (GetRegister(rBX) <= valueToCompare)  getIP().Advance();
  
  return true;
}

bool cHardwareBCR::Inst_IfEquX(cAvidaContext&)       // Execute next if BX == X; X value set according to NOP label
{
  // Compares value in BX to a specific value.  The value to compare to is determined by the nop label as follows:
  //    no nop label (default): valueToCompare = 1;
  //    nop-A: toggles valueToCompare sign-bit 
  //    nop-B: valueToCompare left-shift by 1-bit
  //    nop-C: valueToCompare left-shift by 2-bits
  //    nop-D: valueToCompare left-shift by 3-bits, etc.
  
  int valueToCompare = 1;
  
  ReadLabel();
  const cCodeLabel& shift_label = GetLabel();
  for (int i = 0; i < shift_label.GetSize(); i++) {
    if (shift_label[i] == rAX) {
      valueToCompare *= -1;
    } else {
      valueToCompare <<= shift_label[i];
    }
  }
  
  if (GetRegister(rBX) != valueToCompare)  getIP().Advance();
  
  return true;
}

bool cHardwareBCR::Inst_Pop(cAvidaContext&)
{
  const int reg_used = FindModifiedRegister(rBX);
  DataValue pop = stackPop();
  setInternalValue(reg_used, pop.value, pop);
  return true;
}

bool cHardwareBCR::Inst_Push(cAvidaContext&)
{
  const int reg_used = FindModifiedRegister(rBX);
  getStack(m_threads[m_cur_thread].cur_stack).Push(m_threads[m_cur_thread].reg[reg_used]);
  return true;
}

bool cHardwareBCR::Inst_PopAll(cAvidaContext&)
{
  int reg_used = FindModifiedRegister(rBX);
  for (int i = 0; i < NUM_REGISTERS; i++) {
    DataValue pop = stackPop();
    setInternalValue(reg_used, pop.value, pop);
    reg_used++;
    if (reg_used == NUM_REGISTERS) reg_used = 0;
  }
  return true;
}

bool cHardwareBCR::Inst_PushAll(cAvidaContext&)
{
  int reg_used = FindModifiedRegister(rBX);
  for (int i = 0; i < NUM_REGISTERS; i++) {
    getStack(m_threads[m_cur_thread].cur_stack).Push(m_threads[m_cur_thread].reg[reg_used]);
    reg_used++;
    if (reg_used == NUM_REGISTERS) reg_used = 0;
  }
  return true;
}

bool cHardwareBCR::Inst_SwitchStack(cAvidaContext&) { switchStack(); return true; }

bool cHardwareBCR::Inst_SwapStackTop(cAvidaContext&)
{
  DataValue v0 = getStack(0).Pop();
  DataValue v1 = getStack(1).Pop();
  getStack(0).Push(v1);
  getStack(1).Push(v0);
  return true;
}

bool cHardwareBCR::Inst_Swap(cAvidaContext&)
{
  const int op1 = FindModifiedRegister(rBX);
  const int op2 = FindModifiedNextRegister(op1);
  DataValue v1 = m_threads[m_cur_thread].reg[op1];
  m_threads[m_cur_thread].reg[op1] = m_threads[m_cur_thread].reg[op2];
  m_threads[m_cur_thread].reg[op2] = v1;
  return true;
}

bool cHardwareBCR::Inst_ShiftR(cAvidaContext&)
{
  const int reg_used = FindModifiedRegister(rBX);
  setInternalValue(reg_used, m_threads[m_cur_thread].reg[reg_used].value >> 1,
      m_threads[m_cur_thread].reg[reg_used]);
  return true;
}

bool cHardwareBCR::Inst_ShiftL(cAvidaContext&)
{
  const int reg_used = FindModifiedRegister(rBX);
  setInternalValue(reg_used, m_threads[m_cur_thread].reg[reg_used].value << 1,
      m_threads[m_cur_thread].reg[reg_used]);
  return true;
}


bool cHardwareBCR::Inst_Inc(cAvidaContext&)
{
  const int reg_used = FindModifiedRegister(rBX);
  setInternalValue(reg_used, m_threads[m_cur_thread].reg[reg_used].value + 1,
      m_threads[m_cur_thread].reg[reg_used]);
  return true;
}

bool cHardwareBCR::Inst_Dec(cAvidaContext&)
{
  const int reg_used = FindModifiedRegister(rBX);
  setInternalValue(reg_used, m_threads[m_cur_thread].reg[reg_used].value - 1,
      m_threads[m_cur_thread].reg[reg_used]);
  return true;
}

bool cHardwareBCR::Inst_Zero(cAvidaContext&)
{
  const int reg_used = FindModifiedRegister(rBX);
  setInternalValue(reg_used, 0, false);
  return true;
}

bool cHardwareBCR::Inst_One(cAvidaContext&)
{
  const int reg_used = FindModifiedRegister(rBX);
  setInternalValue(reg_used, 1, false);
  return true;
}

bool cHardwareBCR::Inst_Rand(cAvidaContext&)
{
  const int reg_used = FindModifiedRegister(rBX);
  int randsign = m_world->GetRandom().GetUInt(0,2) ? -1 : 1;
  setInternalValue(reg_used, m_world->GetRandom().GetInt(INT_MAX) * randsign, false);
  return true;
}

bool cHardwareBCR::Inst_Add(cAvidaContext&)
{
  const int dst = FindModifiedRegister(rBX);
  const int op1 = FindModifiedRegister(dst);
  const int op2 = FindModifiedNextRegister(op1);
  DataValue& r1 = m_threads[m_cur_thread].reg[op1];
  DataValue& r2 = m_threads[m_cur_thread].reg[op2];
  setInternalValue(dst, r1.value + r2.value, r1, r2);
  return true;
}

bool cHardwareBCR::Inst_Sub(cAvidaContext&)
{
  const int dst = FindModifiedRegister(rBX);
  const int op1 = FindModifiedRegister(dst);
  const int op2 = FindModifiedNextRegister(op1);
  DataValue& r1 = m_threads[m_cur_thread].reg[op1];
  DataValue& r2 = m_threads[m_cur_thread].reg[op2];
  setInternalValue(dst, r1.value - r2.value, r1, r2);
  return true;
}

bool cHardwareBCR::Inst_Mult(cAvidaContext&)
{
  const int dst = FindModifiedRegister(rBX);
  const int op1 = FindModifiedRegister(dst);
  const int op2 = FindModifiedNextRegister(op1);
  DataValue& r1 = m_threads[m_cur_thread].reg[op1];
  DataValue& r2 = m_threads[m_cur_thread].reg[op2];
  setInternalValue(dst, r1.value * r2.value, r1, r2);
  return true;
}

bool cHardwareBCR::Inst_Div(cAvidaContext&)
{
  const int dst = FindModifiedRegister(rBX);
  const int op1 = FindModifiedRegister(dst);
  const int op2 = FindModifiedNextRegister(op1);
  DataValue& r1 = m_threads[m_cur_thread].reg[op1];
  DataValue& r2 = m_threads[m_cur_thread].reg[op2];
  if (r2.value != 0) {
    if (0 - INT_MAX > r1.value && r2.value == -1)
      m_organism->Fault(FAULT_LOC_MATH, FAULT_TYPE_ERROR, "div: Float exception");
    else
      setInternalValue(dst, r1.value / r2.value, r1, r2);
  } else {
    m_organism->Fault(FAULT_LOC_MATH, FAULT_TYPE_ERROR, "div: dividing by 0");
    return false;
  }
  return true;
}

bool cHardwareBCR::Inst_Mod(cAvidaContext&)
{
  const int dst = FindModifiedRegister(rBX);
  const int op1 = FindModifiedRegister(dst);
  const int op2 = FindModifiedNextRegister(op1);
  DataValue& r1 = m_threads[m_cur_thread].reg[op1];
  DataValue& r2 = m_threads[m_cur_thread].reg[op2];
  if (r2.value != 0) {
    setInternalValue(dst, r1.value % r2.value, r1, r2);
  } else {
    m_organism->Fault(FAULT_LOC_MATH, FAULT_TYPE_ERROR, "mod: modding by 0");
    return false;
  }
  return true;
}


bool cHardwareBCR::Inst_Nand(cAvidaContext&)
{
  const int dst = FindModifiedRegister(rBX);
  const int op1 = FindModifiedRegister(dst);
  const int op2 = FindModifiedNextRegister(op1);
  DataValue& r1 = m_threads[m_cur_thread].reg[op1];
  DataValue& r2 = m_threads[m_cur_thread].reg[op2];
  setInternalValue(dst, ~(r1.value & r2.value), r1, r2);
  return true;
}

bool cHardwareBCR::Inst_SetMemory(cAvidaContext& ctx)
{
  int mem_label = FindModifiedRegister(rBX);
  int head_used = FindModifiedHead(hFLOW);
  
  int mem_id = m_mem_ids[mem_label];
  
  // Check for existing mem_space
  if (m_mem_ids[mem_label] < 0) {
    mem_id = m_mem_array.GetSize();
    m_mem_array.Resize(mem_id + 1);
    m_mem_ids[mem_label] = mem_id;
  }
  
  m_threads[m_cur_thread].heads[head_used].Set(0, mem_id);
  return true;
}

bool cHardwareBCR::Inst_TaskIO(cAvidaContext& ctx)
{
  const int reg_used = FindModifiedRegister(rBX);
  DataValue& reg = m_threads[m_cur_thread].reg[reg_used];
  
  // Do the "put" component
  m_organism->DoOutput(ctx, reg.value);  // Check for tasks completed.
  m_last_output = m_cycle_count;
  
  // Do the "get" component
  const int value_in = m_organism->GetNextInput();
  setInternalValue(reg_used, value_in, true);
  m_organism->DoInput(value_in);
  
  return true;
}

bool cHardwareBCR::Inst_TaskInput(cAvidaContext&)
{
  const int reg_used = FindModifiedRegister(rBX);
  
  // Do the "get" component
  const int value_in = m_organism->GetNextInput();
  setInternalValue(reg_used, value_in, true);
  m_organism->DoInput(value_in);
  
  return true;
}

bool cHardwareBCR::Inst_TaskOutput(cAvidaContext& ctx)
{
  const int reg_used = FindModifiedRegister(rBX);
  DataValue& reg = m_threads[m_cur_thread].reg[reg_used];
  
  // Do the "put" component
  m_organism->DoOutput(ctx, reg.value);  // Check for tasks completed.
  m_last_output = m_cycle_count;
  
  return true;
}

bool cHardwareBCR::Inst_SGMove(cAvidaContext&)
{
  assert(m_ext_mem.GetSize() > 3);
  
  const cStateGrid& sg = m_organism->GetStateGrid();
  
  int& x = m_ext_mem[0];
  int& y = m_ext_mem[1];
  
  const int facing = m_ext_mem[2];
  
  // State grid is treated as a 2-dimensional toroidal grid with size [0, width) and [0, height)
  switch (facing) {
    case 0: // N
      if (++y == sg.GetHeight()) y = 0;
      break;
      
    case 1: // NE
      if (++x == sg.GetWidth()) x = 0;
      if (++y == sg.GetHeight()) y = 0;
      break;
      
    case 2: // E
      if (++x == sg.GetWidth()) x = 0;
      break;
      
    case 3: // SE
      if (++x == sg.GetWidth()) x = 0;
      if (--y == -1) y = sg.GetHeight() - 1;
      break;
      
    case 4: // S
      if (--y == -1) y = sg.GetHeight() - 1;
      break;
      
    case 5: // SW
      if (--x == -1) x = sg.GetWidth() - 1;
      if (--y == -1) y = sg.GetHeight() - 1;
      break;
      
    case 6: // W
      if (--x == -1) x = sg.GetWidth() - 1;
      break;
      
    case 7: // NW
      if (--x == -1) x = sg.GetWidth() - 1;
      if (++y == sg.GetHeight()) y = 0;
      break;
      
    default:
      assert(facing >= 0 && facing <= 7);
  }
  
  // Increment state observed count
  m_ext_mem[3 + sg.GetStateAt(x, y)]++;
  
  // Save this location in the movement history
  m_ext_mem.Push(sg.GetIDFor(x, y));
  return true;
}

bool cHardwareBCR::Inst_SGRotateL(cAvidaContext&)
{
  assert(m_ext_mem.GetSize() > 3);
  if (--m_ext_mem[2] < 0) m_ext_mem[2] = 7;
  return true;
}

bool cHardwareBCR::Inst_SGRotateR(cAvidaContext&)
{
  assert(m_ext_mem.GetSize() > 3);
  if (++m_ext_mem[2] > 7) m_ext_mem[2] = 0;
  return true;
}

bool cHardwareBCR::Inst_SGSense(cAvidaContext&)
{
  const cStateGrid& sg = m_organism->GetStateGrid();
  const int reg_used = FindModifiedRegister(rBX);
  setInternalValue(reg_used, sg.SenseStateAt(m_ext_mem[0], m_ext_mem[1]));
  return true;
}

bool cHardwareBCR::Inst_MoveHead(cAvidaContext&)
{
  const int head_used = FindModifiedHead(hIP);
  const int target = FindModifiedHead(hFLOW);
  getHead(head_used).Set(getHead(target));
  if (head_used == hIP) m_advance_ip = false;
  return true;
}

bool cHardwareBCR::Inst_MoveHeadIfNEqu(cAvidaContext&)
{
  const int op1 = FindModifiedRegister(rBX);
  const int op2 = FindModifiedNextRegister(op1);
  const int head_used = FindModifiedHead(hIP);
  const int target = FindModifiedHead(hFLOW);
  if (m_threads[m_cur_thread].reg[op1].value != m_threads[m_cur_thread].reg[op2].value) {
    getHead(head_used).Set(getHead(target));
    if (head_used == hIP) m_advance_ip = false;
  }
  return true;
}

bool cHardwareBCR::Inst_MoveHeadIfLess(cAvidaContext&)
{
  const int op1 = FindModifiedRegister(rBX);
  const int op2 = FindModifiedNextRegister(op1);
  const int head_used = FindModifiedHead(hIP);
  const int target = FindModifiedHead(hFLOW);
  if (m_threads[m_cur_thread].reg[op1].value < m_threads[m_cur_thread].reg[op2].value) {
    getHead(head_used).Set(getHead(target));
    if (head_used == hIP) m_advance_ip = false;
  }
  return true;
}


bool cHardwareBCR::Inst_JumpHead(cAvidaContext&)
{
  const int head_used = FindModifiedHead(hIP);
  const int reg = FindModifiedRegister(rCX);
  getHead(head_used).Jump(m_threads[m_cur_thread].reg[reg].value);
  if (head_used == hIP) m_advance_ip = false;
  return true;
}

bool cHardwareBCR::Inst_GetHead(cAvidaContext&)
{
  const int head_used = FindModifiedHead(hIP);
  const int reg = FindModifiedRegister(rCX);
  setInternalValue(reg, getHead(head_used).GetPosition());
  return true;
}

bool cHardwareBCR::Inst_IfCopiedCompLabel(cAvidaContext&)
{
  ReadLabel();
  GetLabel().Rotate(1, NUM_NOPS);
  if (GetLabel() != GetReadLabel())  getIP().Advance();
  return true;
}

bool cHardwareBCR::Inst_IfCopiedDirectLabel(cAvidaContext&)
{
  ReadLabel();
  if (GetLabel() != GetReadLabel())  getIP().Advance();
  return true;
}

bool cHardwareBCR::Inst_IfCopiedCompSeq(cAvidaContext&)
{
  ReadLabel();
  GetLabel().Rotate(1, NUM_NOPS);
  if (GetLabel() != GetReadSequence())  getIP().Advance();
  return true;
}

bool cHardwareBCR::Inst_IfCopiedDirectSeq(cAvidaContext&)
{
  ReadLabel();
  if (GetLabel() != GetReadSequence())  getIP().Advance();
  return true;
}


bool cHardwareBCR::Inst_DidCopyCompLabel(cAvidaContext&)
{
  ReadLabel();
  GetLabel().Rotate(1, NUM_NOPS);
  setInternalValue(rBX, (GetLabel() == GetReadLabel()), false);
  return true;
}

bool cHardwareBCR::Inst_DidCopyDirectLabel(cAvidaContext&)
{
  ReadLabel();
  setInternalValue(rBX, (GetLabel() == GetReadLabel()), false);
  return true;
}

bool cHardwareBCR::Inst_DidCopyCompSeq(cAvidaContext&)
{
  ReadLabel();
  GetLabel().Rotate(1, NUM_NOPS);
  setInternalValue(rBX, (GetLabel() == GetReadSequence()), false);
  return true;
}

bool cHardwareBCR::Inst_DidCopyDirectSeq(cAvidaContext&)
{
  ReadLabel();
  setInternalValue(rBX, (GetLabel() == GetReadSequence()), false);
  return true;
}


bool cHardwareBCR::Inst_Divide(cAvidaContext& ctx)
{
  const int mem_space_used = GetHead(hWRITE).GetMemSpace();
  const int write_head_pos = GetHead(hWRITE).GetPosition();

  return Divide_Main(ctx, mem_space_used, write_head_pos, 1.0);
}

bool cHardwareBCR::Inst_DivideMemory(cAvidaContext& ctx)
{
  int mem_space_used = FindModifiedRegister(rBX);

  if (mem_space_used < rBX || m_mem_ids[mem_space_used]  < 0) return false;
  
  mem_space_used = m_mem_ids[mem_space_used];
  int end_of_memory = m_mem_array[mem_space_used].GetSize() - 1;

  return Divide_Main(ctx, mem_space_used, end_of_memory, 1.0);
}

bool cHardwareBCR::Inst_HeadRead(cAvidaContext& ctx)
{
  const int head_id = FindModifiedHead(hREAD);
  const int dst = FindModifiedRegister(rAX);  
  getHead(head_id).Adjust();
  
  // Mutations only occur on the read, for the moment.
  Instruction read_inst;
  if (m_organism->TestCopyMut(ctx)) {
    read_inst = m_inst_set->GetRandomInst(ctx);
  } else {
    read_inst = getHead(head_id).GetInst();
  }
  setInternalValue(dst, read_inst.GetOp());
  ReadInst(read_inst);
  
  if (m_slip_read_head && m_organism->TestCopySlip(ctx))
    getHead(head_id).Set(ctx.GetRandom().GetInt(getHead(head_id).GetMemory().GetSize()));
  
  getHead(head_id).Advance();
  return true;
}

bool cHardwareBCR::Inst_HeadWrite(cAvidaContext& ctx)
{
  const int head_id = FindModifiedHead(hWRITE);
  const int src = FindModifiedRegister(rAX);
  cHeadCPU& active_head = getHead(head_id);
  
  int mem_space_used = active_head.GetMemSpace();
  
  if (active_head.GetPosition() >= m_mem_array[mem_space_used].GetSize() - 1) {
		m_mem_array[mem_space_used].Resize(m_mem_array[mem_space_used].GetSize() + 1);
		m_mem_array[mem_space_used].Copy(m_mem_array[mem_space_used].GetSize() - 1, m_mem_array[mem_space_used].GetSize() - 2);
	}
  
  active_head.Adjust();
  
  int value = m_threads[m_cur_thread].reg[src].value;
  if (value < 0 || value >= m_inst_set->GetSize()) value = 0;
  
  active_head.SetInst(Instruction(value));
  active_head.SetFlagCopied();
  
  if (m_organism->TestCopyIns(ctx)) active_head.InsertInst(m_inst_set->GetRandomInst(ctx));
  if (m_organism->TestCopyDel(ctx)) active_head.RemoveInst();
  if (m_organism->TestCopyUniform(ctx)) doUniformCopyMutation(ctx, active_head);
  if (!m_slip_read_head && m_organism->TestCopySlip(ctx)) 
    doSlipMutation(ctx, active_head.GetMemory(), active_head.GetPosition());
  
  // Advance the head after write...
  active_head++;
  return true;
}

bool cHardwareBCR::Inst_HeadCopy(cAvidaContext& ctx)
{
  // For the moment, this cannot be nop-modified.
  cHeadCPU& read_head = getHead(hREAD);
  cHeadCPU& write_head = getHead(hWRITE);
  
  int mem_space_used = write_head.GetMemSpace();
  
  if (write_head.GetPosition() >= m_mem_array[mem_space_used].GetSize() - 1) {
		m_mem_array[mem_space_used].Resize(m_mem_array[mem_space_used].GetSize() + 1);
		m_mem_array[mem_space_used].Copy(m_mem_array[mem_space_used].GetSize() - 1, m_mem_array[mem_space_used].GetSize() - 2);
	}

  
  read_head.Adjust();
  write_head.Adjust();
  
  // Do mutations.
  Instruction read_inst = read_head.GetInst();
  ReadInst(read_inst);
  if (m_organism->TestCopyMut(ctx)) {
    read_inst = m_inst_set->GetRandomInst(ctx);
    write_head.SetFlagMutated();
    write_head.SetFlagCopyMut();
  }
  
  write_head.SetInst(read_inst);
  write_head.SetFlagCopied();  // Set the copied flag...
  
  if (m_organism->TestCopyIns(ctx)) write_head.InsertInst(m_inst_set->GetRandomInst(ctx));
  if (m_organism->TestCopyDel(ctx)) write_head.RemoveInst();
  if (m_organism->TestCopyUniform(ctx)) doUniformCopyMutation(ctx, write_head);
  if (m_organism->TestCopySlip(ctx)) {
    if (m_slip_read_head) {
      read_head.Set(ctx.GetRandom().GetInt(read_head.GetMemory().GetSize()));
    } else 
      doSlipMutation(ctx, write_head.GetMemory(), write_head.GetPosition());
  }
  
  read_head.Advance();
  write_head.Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Label_Comp_S(cAvidaContext&)
{
  ReadLabel();
  GetLabel().Rotate(1, NUM_NOPS);
  FindLabelStart(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Label_Comp_F(cAvidaContext&)
{
  ReadLabel();
  GetLabel().Rotate(1, NUM_NOPS);
  FindLabelForward(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Label_Comp_B(cAvidaContext&)
{
  ReadLabel();
  GetLabel().Rotate(1, NUM_NOPS);
  FindLabelBackward(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Label_Direct_S(cAvidaContext&)
{
  ReadLabel();
  FindLabelStart(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Label_Direct_F(cAvidaContext&)
{
  ReadLabel();
  FindLabelForward(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Label_Direct_B(cAvidaContext&)
{
  ReadLabel();
  FindLabelBackward(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Seq_Comp_S(cAvidaContext&)
{
  ReadLabel();
  GetLabel().Rotate(1, NUM_NOPS);
  FindNopSequenceStart(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Seq_Comp_F(cAvidaContext&)
{
  ReadLabel();
  GetLabel().Rotate(1, NUM_NOPS);
  FindNopSequenceForward(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Seq_Comp_B(cAvidaContext&)
{
  ReadLabel();
  GetLabel().Rotate(1, NUM_NOPS);
  FindNopSequenceBackward(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Seq_Direct_S(cAvidaContext&)
{
  ReadLabel();
  FindNopSequenceStart(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Seq_Direct_F(cAvidaContext&)
{
  ReadLabel();
  FindNopSequenceForward(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_Search_Seq_Direct_B(cAvidaContext&)
{
  ReadLabel();
  FindNopSequenceBackward(getHead(hFLOW), true);
  getHead(hFLOW).Advance();
  return true;
}

bool cHardwareBCR::Inst_WaitCondition_Equal(cAvidaContext&)
{
  const int wait_value = FindModifiedRegister(rBX);
  const int check_reg = FindModifiedRegister(rDX);
  const int wait_dst = FindModifiedRegister(wait_value);
  
  // Check if condition has already been met
  for (int i = 0; i < m_threads.GetSize(); i++) {
    if (i != m_cur_thread && m_threads[i].reg[check_reg].value == m_threads[m_cur_thread].reg[wait_value].value) {
      setInternalValue(wait_dst, m_threads[i].reg[check_reg].value, m_threads[i].reg[check_reg]);
      return true;
    }
  }
  
  // Fail to sleep if this is the last thread awake
  if (m_waiting_threads == m_running_threads) return false;
  
  // Put thread to sleep with appropriate wait condition
  m_threads[m_cur_thread].active = false;
  m_waiting_threads++;
  m_threads[m_cur_thread].wait_equal = true;
  m_threads[m_cur_thread].wait_less = false;
  m_threads[m_cur_thread].wait_greater = false;
  m_threads[m_cur_thread].wait_reg = check_reg;
  m_threads[m_cur_thread].wait_value = m_threads[m_cur_thread].reg[wait_value].value;
  m_threads[m_cur_thread].wait_dst = wait_dst;
  
  return true;
}

bool cHardwareBCR::Inst_WaitCondition_Less(cAvidaContext&)
{
  const int wait_value = FindModifiedRegister(rBX);
  const int check_reg = FindModifiedRegister(rDX);
  const int wait_dst = FindModifiedRegister(wait_value);
  
  // Check if condition has already been met
  for (int i = 0; i < m_threads.GetSize(); i++) {
    if (i != m_cur_thread && m_threads[i].reg[check_reg].value <
            m_threads[m_cur_thread].reg[wait_value].value) {
      setInternalValue(wait_dst, m_threads[i].reg[check_reg].value,
                m_threads[i].reg[check_reg]);
      return true;
    }
  }
  
  // Fail to sleep if this is the last thread awake
  if (m_waiting_threads == m_running_threads) return false;
  
  // Put thread to sleep with appropriate wait condition
  m_threads[m_cur_thread].active = false;
  m_waiting_threads++;
  m_threads[m_cur_thread].wait_equal = false;
  m_threads[m_cur_thread].wait_less = true;
  m_threads[m_cur_thread].wait_greater = false;
  m_threads[m_cur_thread].wait_reg = check_reg;
  m_threads[m_cur_thread].wait_value = m_threads[m_cur_thread].reg[wait_value].value;
  m_threads[m_cur_thread].wait_dst = wait_dst;
  
  return true;
}

bool cHardwareBCR::Inst_WaitCondition_Greater(cAvidaContext&)
{
  const int wait_value = FindModifiedRegister(rBX);
  const int check_reg = FindModifiedRegister(rDX);
  const int wait_dst = FindModifiedRegister(wait_value);
  
  // Check if condition has already been met
  for (int i = 0; i < m_threads.GetSize(); i++) {
    if (i != m_cur_thread && m_threads[i].reg[check_reg].value >
            m_threads[m_cur_thread].reg[wait_value].value) {
      setInternalValue(wait_dst, m_threads[i].reg[check_reg].value,
              m_threads[i].reg[check_reg]);
      return true;
    }
  }
  
  // Fail to sleep if this is the last thread awake
  if (m_waiting_threads == m_running_threads) return false;
  
  // Put thread to sleep with appropriate wait condition
  m_threads[m_cur_thread].active = false;
  m_waiting_threads++;
  m_threads[m_cur_thread].wait_equal = false;
  m_threads[m_cur_thread].wait_less = false;
  m_threads[m_cur_thread].wait_greater = true;
  m_threads[m_cur_thread].wait_reg = check_reg;
  m_threads[m_cur_thread].wait_value = m_threads[m_cur_thread].reg[wait_value].value;
  m_threads[m_cur_thread].wait_dst = wait_dst;
  
  return true;
}

bool cHardwareBCR::Inst_Repro(cAvidaContext& ctx)
{
  // these checks should be done, but currently they make some assumptions
  // that crash when evaluating this kind of organism -- JEB
  
  cCPUMemory& memory = m_mem_array[0];
  
  if (m_organism->GetPhenotype().GetCurBonus() < m_world->GetConfig().REQUIRED_BONUS.Get()) return false;
  
  // Since the divide will now succeed, set up the information to be sent
  // to the new organism
  InstructionSequencePtr offspring_seq(new InstructionSequence(memory));
  HashPropertyMap props;
  cHardwareManager::SetupPropertyMap(props, (const char*)m_inst_set->GetInstSetName());
  Genome offspring(GetType(), props, offspring_seq);

  m_organism->OffspringGenome() = offspring;  
  m_organism->GetPhenotype().SetLinesCopied(memory.GetSize());
  
  int lines_executed = 0;
  for (int i = 0; i < memory.GetSize(); i++) if (memory.FlagExecuted(i)) lines_executed++;
  m_organism->GetPhenotype().SetLinesExecuted(lines_executed);
  
  const Genome& org = m_organism->GetGenome();
  ConstInstructionSequencePtr org_seq_p;
  org_seq_p.DynamicCastFrom(org.Representation());
  const InstructionSequence& org_genome = *org_seq_p;  

  Genome& child = m_organism->OffspringGenome();
  InstructionSequencePtr child_seq_p;
  child_seq_p.DynamicCastFrom(child.Representation());
  InstructionSequence& child_seq = *child_seq_p;  

  // Perform Copy Mutations...
  if (m_organism->GetCopyMutProb() > 0) { // Skip this if no mutations....
    for (int i = 0; i < child_seq.GetSize(); i++) {
//    for (int i = 0; i < m_memory.GetSize(); i++) {
      if (m_organism->TestCopyMut(ctx)) child_seq[i] = m_inst_set->GetRandomInst(ctx);
    }
  }
  
  // Handle Divide Mutations...
  Divide_DoMutations(ctx);
  
  const bool viable = Divide_CheckViable(ctx, org_genome.GetSize(), child_seq.GetSize(), 1);
  if (viable == false) return false;
  
  // Many tests will require us to run the offspring through a test CPU;
  // this is, for example, to see if mutations need to be reverted or if
  // lineages need to be updated.
  Divide_TestFitnessMeasures(ctx);
  
  // reset first time instruction costs
  for (int i = 0; i < m_inst_ft_cost.GetSize(); i++) {
    m_inst_ft_cost[i] = m_inst_set->GetFTCost(Instruction(i));
  }
  
  if (m_world->GetConfig().DIVIDE_METHOD.Get() == DIVIDE_METHOD_SPLIT) {
    m_advance_ip = false;
  }
  
  // Activate the child
  bool parent_alive = m_organism->ActivateDivide(ctx);
  
  // Do more work if the parent lives through the birth of the offspring
  if (parent_alive) {
    if (m_world->GetConfig().DIVIDE_METHOD.Get() == DIVIDE_METHOD_SPLIT) Reset(ctx);
  }
  
  return true;
}

bool cHardwareBCR::Inst_Die(cAvidaContext& ctx)
{
  m_organism->Die(ctx);
  
  return true;
}

bool cHardwareBCR::Inst_Move(cAvidaContext& ctx)
{
  // In TestCPU, movement fails...
  if (m_organism->GetOrgInterface().GetCellID() == -1) return false;
  
  bool move_success = false;
  if (!m_use_avatar) move_success = m_organism->Move(ctx);
  else if (m_use_avatar) move_success = m_organism->MoveAV(ctx);
  const int out_reg = FindModifiedRegister(rBX);   
  setInternalValue(out_reg, move_success, true);   
  return true;
}

bool cHardwareBCR::Inst_GetNorthOffset(cAvidaContext& ctx) {
  const int out_reg = FindModifiedRegister(rBX);
  int compass_dir = m_organism->GetOrgInterface().GetFacedDir();
  if (m_use_avatar) compass_dir = m_organism->GetOrgInterface().GetAVFacing();
  setInternalValue(out_reg, compass_dir, true);
  return true;
}

bool cHardwareBCR::Inst_GetPositionOffset(cAvidaContext&) {
  const int out_reg = FindModifiedRegister(rBX);
  setInternalValue(out_reg, m_organism->GetNortherly(), true);
  setInternalValue(FindModifiedNextRegister(out_reg), m_organism->GetEasterly(), true);
  return true;
}

bool cHardwareBCR::Inst_GetNortherly(cAvidaContext&) {
  const int out_reg = FindModifiedRegister(rBX);
  setInternalValue(out_reg, m_organism->GetNortherly(), true);
  return true;  
}

bool cHardwareBCR::Inst_GetEasterly(cAvidaContext&) {
  const int out_reg = FindModifiedRegister(rBX);
  setInternalValue(out_reg, m_organism->GetEasterly(), true);
  return true;  
}

bool cHardwareBCR::Inst_ZeroEasterly(cAvidaContext&) {
  m_organism->ClearEasterly();
  return true;
}

bool cHardwareBCR::Inst_ZeroNortherly(cAvidaContext&) {
  m_organism->ClearNortherly();
  return true;
}

bool cHardwareBCR::Inst_ZeroPosOffset(cAvidaContext&) {
  const int offset = GetRegister(FindModifiedRegister(rBX)) % 3;
  if (offset == 0) {
    m_organism->ClearEasterly();
    m_organism->ClearNortherly();    
  }
  else if (offset == 1) m_organism->ClearEasterly();
  else if (offset == 2) m_organism->ClearNortherly();
  return true;
}

bool cHardwareBCR::Inst_RotateHome(cAvidaContext&)
{
  // Will rotate organism to face birth cell if org never used zero-easterly or zero-northerly. Otherwise will rotate org
  // to face the 'marked' spot where those instructions were executed.
  int easterly = m_organism->GetEasterly();
  int northerly = m_organism->GetNortherly();
  int correct_facing = 0;
  if (northerly > 0 && easterly == 0) correct_facing = 0; // rotate N    
  else if (northerly > 0 && easterly < 0) correct_facing = 1; // rotate NE
  else if (northerly == 0 && easterly < 0) correct_facing = 2; // rotate E
  else if (northerly < 0 && easterly < 0) correct_facing = 3; // rotate SE
  else if (northerly < 0 && easterly == 0) correct_facing = 4; // rotate S
  else if (northerly < 0 && easterly > 0) correct_facing = 5; // rotate SW
  else if (northerly == 0 && easterly > 0) correct_facing = 6; // rotate W
  else if (northerly > 0 && easterly > 0) correct_facing = 7; // rotate NW  
  
  int rotates = m_organism->GetNeighborhoodSize();
  if (m_use_avatar == 2) rotates = m_organism->GetOrgInterface().GetAVNumNeighbors();
  for (int i = 0; i < rotates; i++) {
    m_organism->Rotate(1);
    if (!m_use_avatar && m_organism->GetOrgInterface().GetFacedDir() == correct_facing) break;
    else if (m_use_avatar && m_organism->GetOrgInterface().GetAVFacing() == correct_facing) break;
  }
  return true;
}

bool cHardwareBCR::Inst_RotateUnoccupiedCell(cAvidaContext&)
{
  if (m_use_avatar && m_use_avatar != 2) return false;
  const int reg_used = FindModifiedRegister(rBX);
  
  int num_neighbors = m_organism->GetNeighborhoodSize();
  if (m_use_avatar) num_neighbors = m_organism->GetOrgInterface().GetAVNumNeighbors();
  for (int i = 0; i < num_neighbors; i++) {
    if ((!m_use_avatar && !m_organism->IsNeighborCellOccupied()) || (m_use_avatar == 2 && !m_organism->GetOrgInterface().FacedHasAV())) { 
      setInternalValue(reg_used, 1, true);      
      return true;
    }
    m_organism->Rotate(1); // continue to rotate
  }  
  setInternalValue(reg_used, 0, true);
  return true;
}

bool cHardwareBCR::Inst_RotateX(cAvidaContext&)
{
  int num_neighbors = m_organism->GetNeighborhoodSize();
  if (m_use_avatar) num_neighbors = m_organism->GetOrgInterface().GetAVNumNeighbors();
  int rot_dir = 1;
  // If this organism has no neighbors, ignore rotate.
  if (num_neighbors == 0) return false;
  
  const int reg_used = FindModifiedRegister(rBX);
  int rot_num = m_threads[m_cur_thread].reg[reg_used].value;
  // rotate the nop number of times in the appropriate direction
  rot_num < 0 ? rot_dir = -1 : rot_dir = 1;
  rot_num = abs(rot_num);
  if (rot_num > 7) rot_num = rot_num % 8;
  for (int i = 0; i < rot_num; i++) m_organism->Rotate(rot_dir);
  setInternalValue(reg_used, rot_num * rot_dir, true);
  return true;
}

bool cHardwareBCR::Inst_SenseResourceID(cAvidaContext& ctx)
{
  Apto::Array<double> cell_res;
  if (!m_use_avatar) cell_res = m_organism->GetOrgInterface().GetResources(ctx);
  else if (m_use_avatar) cell_res = m_organism->GetOrgInterface().GetAVResources(ctx); 
  int reg_to_set = FindModifiedRegister(rBX);  
  double max_resource = 0.0;    
  // if more than one resource is available, return the resource ID with the most available in this spot (note that, with global resources, the GLOBAL total will evaluated)
  for (int i = 0; i < cell_res.GetSize(); i++) {
    if (cell_res[i] > max_resource) {
      max_resource = cell_res[i];
      setInternalValue(reg_to_set, i, true);
    }
  }    
  return true;
}

bool cHardwareBCR::Inst_SenseNest(cAvidaContext& ctx)
{
  Apto::Array<double> cell_res;
  if (!m_use_avatar) cell_res = m_organism->GetOrgInterface().GetResources(ctx);
  else if (m_use_avatar) cell_res = m_organism->GetOrgInterface().GetAVResources(ctx); 
  
  const cResourceLib& resource_lib = m_world->GetEnvironment().GetResourceLib();
  const int reg_used = FindModifiedRegister(rBX);
  
  int nest_id = m_threads[m_cur_thread].reg[reg_used].value;
  int nest_val = 0;
  
  // if invalid nop value, return the id of the first nest in the cell with val >= 1
  if (nest_id < 0 || nest_id >= resource_lib.GetSize() || resource_lib.GetResource(nest_id)->GetHabitat() != 3) {
    for (int i = 0; i < cell_res.GetSize(); i++) {
      if (resource_lib.GetResource(i)->GetHabitat() == 3 && cell_res[i] >= 1) {
        nest_id = i;
        nest_val = (int) cell_res[i];
        break;
      }
    }
  }
  else nest_val = (int) cell_res[nest_id];  
  
  setInternalValue(reg_used, nest_id, true);
  const int val_reg = FindModifiedNextRegister(reg_used);
  setInternalValue(val_reg, nest_val, true);
  return true;
}

bool cHardwareBCR::Inst_LookAhead(cAvidaContext& ctx)
{
  int cell = m_organism->GetOrgInterface().GetCellID();
  int facing = m_organism->GetOrgInterface().GetFacedDir();
  if (m_use_avatar) {
    cell = m_organism->GetOrgInterface().GetAVCellID();
    facing = m_organism->GetOrgInterface().GetAVFacing();
  }
  return GoLook(ctx, facing, cell);
}

// Will return relative org facing (rotations to intercept) rather than group info for sighted org
bool cHardwareBCR::Inst_LookAheadIntercept(cAvidaContext& ctx)
{
  m_sensor.SetReturnRelativeFacing(true);
  return Inst_LookAhead(ctx);
}

bool cHardwareBCR::Inst_LookAround(cAvidaContext& ctx)
{
  // dir register is 5th mod (will be count reg)
  int hab_reg = FindModifiedRegister(rBX);
  int dist_reg = FindModifiedNextRegister(hab_reg);
  int st_reg = FindModifiedNextRegister(dist_reg);
  int id_reg = FindModifiedNextRegister(st_reg);
  int dir_reg = FindModifiedNextRegister(id_reg);
  
  int search_dir = abs(m_threads[m_cur_thread].reg[dir_reg].value) % 3;
  
  if (m_world->GetConfig().LOOK_DISABLE.Get() == 5) {
    int org_type = m_world->GetConfig().LOOK_DISABLE_TYPE.Get();
    bool is_target_type = false;
    if (org_type == 0 && m_organism->GetForageTarget() <= -2) is_target_type = true;
    else if (org_type == 1 && m_organism->GetForageTarget() >-2) is_target_type = true;
    else if (org_type == 2) is_target_type = true;
    
    if (is_target_type) {
      int rand = m_world->GetRandom().GetInt(INT_MAX);
      search_dir = rand % 3;
    }
  }
  
  if (search_dir == 1) search_dir = -1;
  else if (search_dir == 2) search_dir = 1;
  
  int facing = m_organism->GetOrgInterface().GetFacedDir() + search_dir;
  if (m_use_avatar) facing = m_organism->GetOrgInterface().GetAVFacing() + search_dir;
  if (facing == -1) facing = 7;
  else if (facing == 9) facing = 1;
  else if (facing == 8) facing = 0;

  int cell = m_organism->GetOrgInterface().GetCellID();
  if (m_use_avatar) cell = m_organism->GetOrgInterface().GetAVCellID();
  return GoLook(ctx, facing, cell);
}

bool cHardwareBCR::Inst_LookAroundIntercept(cAvidaContext& ctx)
{
  m_sensor.SetReturnRelativeFacing(true);
  return Inst_LookAround(ctx);
}

bool cHardwareBCR::Inst_LookFT(cAvidaContext& ctx)
{
  // override any org inputs and just let this org see the food resource that matches it's forage target (not designed for predators)
  int cell = m_organism->GetOrgInterface().GetCellID();
  int facing = m_organism->GetOrgInterface().GetFacedDir();
  if (m_use_avatar) { 
    facing = m_organism->GetOrgInterface().GetAVFacing();
    cell = m_organism->GetOrgInterface().GetAVCellID();
  }
  return GoLook(ctx, facing, cell, true);
}

bool cHardwareBCR::Inst_LookAroundFT(cAvidaContext& ctx)
{
  // dir register is 5th mod (will be count reg)
  int hab_reg = FindModifiedRegister(rBX);
  int dist_reg = FindModifiedNextRegister(hab_reg);
  int st_reg = FindModifiedNextRegister(dist_reg);
  int id_reg = FindModifiedNextRegister(st_reg);
  int dir_reg = FindModifiedNextRegister(id_reg);
  
  int search_dir = abs(m_threads[m_cur_thread].reg[dir_reg].value) % 3;
  
  if (m_world->GetConfig().LOOK_DISABLE.Get() == 5) {
    int org_type = m_world->GetConfig().LOOK_DISABLE_TYPE.Get();
    bool is_target_type = false;
    if (org_type == 0 && m_organism->GetForageTarget() <= -2) is_target_type = true;
    else if (org_type == 1 && m_organism->GetForageTarget() > -2) is_target_type = true;
    else if (org_type == 2) is_target_type = true;
    
    if (is_target_type) {
      int rand = m_world->GetRandom().GetInt(INT_MAX);
      search_dir = rand % 3;
    }
  }
  
  if (search_dir == 1) search_dir = -1;
  else if (search_dir == 2) search_dir = 1;
  
  int facing = m_organism->GetOrgInterface().GetFacedDir() + search_dir;
  if (m_use_avatar) facing = m_organism->GetOrgInterface().GetAVFacing() + search_dir;
  if (facing == -1) facing = 7;
  else if (facing == 9) facing = 1;
  else if (facing == 8) facing = 0;
  
  int cell = m_organism->GetOrgInterface().GetCellID();
  if (m_use_avatar) cell = m_organism->GetOrgInterface().GetAVCellID();
  return GoLook(ctx, facing, cell, true);
}

bool cHardwareBCR::GoLook(cAvidaContext& ctx, const int look_dir, const int cell_id, bool use_ft) 
{
  // temp check on world geometry until code can handle other geometries
  if (m_world->GetConfig().WORLD_GEOMETRY.Get() != 1) {
    // Instruction sense-diff-ahead only written to work in bounded grids
    return false;
  }  
  if (NUM_REGISTERS < 8) m_world->GetDriver().Feedback().Error("Instruction look-ahead requires at least 8 registers");
  if (!m_use_avatar && m_organism->GetNeighborhoodSize() == 0) return false;
  else if (m_use_avatar && m_organism->GetOrgInterface().GetAVNumNeighbors() == 0) return false;
  
  // define our input (4) and output registers (8)
  sLookRegAssign reg_defs;
  reg_defs.habitat = FindModifiedRegister(rBX);
  // fail if the org is trying to sense a nest/hidden habitat
  int habitat_used = m_threads[m_cur_thread].reg[reg_defs.habitat].value;
  if (habitat_used == 3) return false;
  reg_defs.distance = FindModifiedNextRegister(reg_defs.habitat);
  reg_defs.search_type = FindModifiedNextRegister(reg_defs.distance);
  reg_defs.id_sought = FindModifiedNextRegister(reg_defs.search_type);
  reg_defs.count = FindModifiedNextRegister(reg_defs.id_sought);
  reg_defs.value = FindModifiedNextRegister(reg_defs.count);
  reg_defs.group = FindModifiedNextRegister(reg_defs.value);
  reg_defs.ft = FindModifiedNextRegister(reg_defs.group);
  
  cOrgSensor::sLookOut look_results;
  look_results.report_type = 0;
  look_results.habitat = 0;
  look_results.distance = -1;
  look_results.search_type = 0;
  look_results.id_sought = -1;
  look_results.count = 0;
  look_results.value = 0;
  look_results.group = -9;
  look_results.forage = -9;

  look_results = InitLooking(ctx, reg_defs, look_dir, cell_id, use_ft);
  LookResults(reg_defs, look_results);
  return true;
}

cOrgSensor::sLookOut cHardwareBCR::InitLooking(cAvidaContext& ctx, sLookRegAssign& in_defs, int facing, int cell_id, bool use_ft)
{
  const int habitat_reg = in_defs.habitat;
  const int distance_reg = in_defs.distance;
  const int search_reg = in_defs.search_type;
  const int id_reg = in_defs.id_sought;
  
  cOrgSensor::sLookInit reg_init;
  reg_init.habitat = m_threads[m_cur_thread].reg[habitat_reg].value;
  reg_init.distance = m_threads[m_cur_thread].reg[distance_reg].value;
  reg_init.search_type = m_threads[m_cur_thread].reg[search_reg].value;
  reg_init.id_sought = m_threads[m_cur_thread].reg[id_reg].value;

  return m_sensor.SetLooking(ctx, reg_init, facing, cell_id, use_ft);
}    

void cHardwareBCR::LookResults(sLookRegAssign& regs, cOrgSensor::sLookOut& results)
{
  // habitat_reg=0, distance_reg=1, search_type_reg=2, id_sought_reg=3, count_reg=4, value_reg=5, group_reg=6, forager_type_reg=7
  // return defaults for failed to find
  if (results.report_type == 0) {
    setInternalValue(regs.habitat, results.habitat, true);
    setInternalValue(regs.distance, -1, true);
    setInternalValue(regs.search_type, results.search_type, true);
    setInternalValue(regs.id_sought, results.id_sought, true);
    setInternalValue(regs.count, 0, true);
    setInternalValue(regs.value, 0, true);
    setInternalValue(regs.group, -9, true);
    setInternalValue(regs.ft, -9, true);  
  }
  // report results as sent
  else if (results.report_type == 1) {
    setInternalValue(regs.habitat, results.habitat, true);
    setInternalValue(regs.distance, results.distance, true);
    setInternalValue(regs.search_type, results.search_type, true);
    setInternalValue(regs.id_sought, results.id_sought, true);
    setInternalValue(regs.count, results.count, true);
    setInternalValue(regs.value, results.value, true);
    setInternalValue(regs.group, results.group, true);
    setInternalValue(regs.ft, results.forage, true);  
  }
  
  if (m_world->GetConfig().LOOK_DISABLE.Get() > 5) {
    int org_type = m_world->GetConfig().LOOK_DISABLE_TYPE.Get();
    bool is_target_type = false;
    if (org_type == 0 && m_organism->GetForageTarget() <= -2) is_target_type = true;
    else if (org_type == 1 && m_organism->GetForageTarget() > -2) is_target_type = true;
    else if (org_type == 2) is_target_type = true;
    
    if (is_target_type) {
      int randsign = m_world->GetRandom().GetUInt(0,2) ? -1 : 1;
      int rand = m_world->GetRandom().GetInt(INT_MAX) * randsign;
      int target_reg = m_world->GetConfig().LOOK_DISABLE.Get();
      
      if (target_reg == 6) setInternalValue(regs.habitat, rand, true);
      else if (target_reg == 7) setInternalValue(regs.distance, rand, true);
      else if (target_reg == 8) setInternalValue(regs.search_type, rand, true);
      else if (target_reg == 9) setInternalValue(regs.id_sought, rand, true);
      else if (target_reg == 10) setInternalValue(regs.count, rand, true);
      else if (target_reg == 11) setInternalValue(regs.value, rand, true);
      else if (target_reg == 12) setInternalValue(regs.group, rand, true);
      else if (target_reg == 13) setInternalValue(regs.ft, rand, true);  
    }
  }
  return;
}

bool cHardwareBCR::Inst_SenseFacedHabitat(cAvidaContext& ctx) 
{
  int reg_to_set = FindModifiedRegister(rBX);
  
  // get the resource library
  const cResourceLib& resource_lib = m_world->GetEnvironment().GetResourceLib();
  
  // get the destination cell resource levels
  Apto::Array<double> cell_res;
  if (!m_use_avatar) cell_res = m_organism->GetOrgInterface().GetResources(ctx);
  else if (m_use_avatar) cell_res = m_organism->GetOrgInterface().GetAVResources(ctx); 
  
  // check for any habitats ahead that affect movement, returning the most 'severe' habitat type
  // simulated predator ahead
  for (int i = 0; i < cell_res.GetSize(); i++) {
    if (resource_lib.GetResource(i)->GetHabitat() == 5 && cell_res[i] > 0) {
      setInternalValue(reg_to_set, 3, true);
      return true;
    }    
  }
  // are there any barrier resources in the faced cell    
  for (int i = 0; i < cell_res.GetSize(); i++) {
    if (resource_lib.GetResource(i)->GetHabitat() == 2 && cell_res[i] > 0) {
      setInternalValue(reg_to_set, 2, true);
      return true;
    }    
  }
  // if no barriers, are there any hills in the faced cell    
  for (int i = 0; i < cell_res.GetSize(); i++) {
    if (resource_lib.GetResource(i)->GetHabitat() == 1 && cell_res[i] > 0) {
      setInternalValue(reg_to_set, 1, true);
      return true;
    }
  }
  // if no barriers or hills, we return a 0 to indicate clear sailing
  setInternalValue(reg_to_set, 0, true);
  return true;
}

bool cHardwareBCR::Inst_SetForageTarget(cAvidaContext&)
{
  assert(m_organism != 0);
  int prop_target = GetRegister(FindModifiedRegister(rBX));
  
  //return false if org setting target to current one (avoid paying costs for not switching)
  const int old_target = m_organism->GetForageTarget();
  if (old_target == prop_target) return false;

  // return false if predator trying to become prey and this has been disallowed
  if (old_target <= -2 && prop_target > -2 && (m_world->GetConfig().PRED_PREY_SWITCH.Get() == 0 || m_world->GetConfig().PRED_PREY_SWITCH.Get() == 2)) return false;
  
  // return false if trying to become predator and there are none in the experiment
  if (prop_target <= -2 && m_world->GetConfig().PRED_PREY_SWITCH.Get() < 0) return false;
  
  // return false if trying to become predator this has been disallowed via setforagetarget
  if (prop_target <= -2 && m_world->GetConfig().PRED_PREY_SWITCH.Get() == 2) return false;
  
  // a little mod help...can't set to -1, that's for juevniles only...so only exception to mod help is -2
  if (!m_world->GetEnvironment().IsTargetID(prop_target) && prop_target != -2 && prop_target != -3) {
    int num_fts = 0;
    std::set<int> fts_avail = m_world->GetEnvironment().GetTargetIDs();
    set <int>::iterator itr;    
    for (itr = fts_avail.begin();itr!=fts_avail.end();itr++) if (*itr != -1 && *itr != -2 && *itr != -3) num_fts++;
    if (m_world->GetEnvironment().IsTargetID(-1) && num_fts == 0) prop_target = -1;
    else {
      // ft's may not be sequentially numbered
      int ft_num = abs(prop_target) % num_fts;
      itr = fts_avail.begin();
      for (int i = 0; i < ft_num; i++) itr++;
      prop_target = *itr;
    }
  }

  // make sure we use a valid (resource) target
  // -2 target means setting to predator
  // if (!m_world->GetEnvironment().IsTargetID(prop_target) && (prop_target != -2)) return false;

  // switching between predator and prey means having to switch avatar list...don't run this for orgs with AVCell == -1 (avatars off or test cpu)
  if (m_use_avatar && (((prop_target == -2 || prop_target == -3) && old_target > -2) || (prop_target > -2 && (old_target == -2 || old_target == -3))) &&
      (m_organism->GetOrgInterface().GetAVCellID() != -1)) {
    m_organism->GetOrgInterface().SwitchPredPrey();
    m_organism->SetForageTarget(prop_target);
  }
  else m_organism->SetForageTarget(prop_target);
    
  // Set the new target and return the value
  m_organism->RecordFTSet();
  setInternalValue(FindModifiedRegister(rBX), prop_target, false);
  return true;
}

bool cHardwareBCR::Inst_SetForageTargetOnce(cAvidaContext& ctx)
{
  assert(m_organism != 0);
  if (m_organism->HasSetFT()) return false;
  else return Inst_SetForageTarget(ctx);
}

bool cHardwareBCR::Inst_SetRandForageTargetOnce(cAvidaContext& ctx)
{
  assert(m_organism != 0);
  int cap = 0;
  if (m_world->GetConfig().POPULATION_CAP.Get()) cap = m_world->GetConfig().POPULATION_CAP.Get();
  else if (m_world->GetConfig().POP_CAP_ELDEST.Get()) cap = m_world->GetConfig().POP_CAP_ELDEST.Get();
  if (cap && (m_organism->GetOrgInterface().GetLiveOrgList().GetSize() >= (((double)(cap)) * 0.5)) && m_world->GetRandom().P(0.5)) {
    if (m_organism->HasSetFT()) return false;
    else {
      int num_fts = 0;
      std::set<int> fts_avail = m_world->GetEnvironment().GetTargetIDs();
      set <int>::iterator itr;
      for (itr = fts_avail.begin();itr!=fts_avail.end();itr++) if (*itr != -1 && *itr != -2 && *itr != -3) num_fts++;
      int prop_target = m_world->GetRandom().GetUInt(num_fts);
      if (m_world->GetEnvironment().IsTargetID(-1) && num_fts == 0) prop_target = -1;
      else {
        // ft's may not be sequentially numbered
        int ft_num = abs(prop_target) % num_fts;
        itr = fts_avail.begin();
        for (int i = 0; i < ft_num; i++) itr++;
        prop_target = *itr;
      }
      // Set the new target and return the value
      m_organism->SetForageTarget(prop_target);
      m_organism->RecordFTSet();
      setInternalValue(FindModifiedRegister(rBX), prop_target, false);
      return true;
    }
  }
  else return Inst_SetForageTargetOnce(ctx);
}

bool cHardwareBCR::Inst_GetForageTarget(cAvidaContext& ctx)
{
  assert(m_organism != 0);
  const int target_reg = FindModifiedRegister(rBX);
  setInternalValue(target_reg, m_organism->GetForageTarget(), false);
  return true;
}

bool cHardwareBCR::Inst_CollectSpecific(cAvidaContext& ctx)
{
  const int resource = m_world->GetConfig().COLLECT_SPECIFIC_RESOURCE.Get();
  double res_before = m_organism->GetRBin(resource);
  bool success = DoActualCollect(ctx, resource, false);
  double res_after = m_organism->GetRBin(resource);
  int out_reg = FindModifiedRegister(rBX);
  setInternalValue(out_reg, (int)(res_after - res_before), true);
  setInternalValue(FindModifiedNextRegister(out_reg), (int)(res_after), true);
  return success;
}

bool cHardwareBCR::Inst_GetResStored(cAvidaContext& ctx)
{
  int resource_id = abs(GetRegister(FindModifiedRegister(rBX)));
  Apto::Array<double> bins = m_organism->GetRBins();
  resource_id %= bins.GetSize();
  int out_reg = FindModifiedRegister(rBX);
  setInternalValue(out_reg, (int)(bins[resource_id]), true);
  return true;
}

// Sets organism's opinion to the value in ?BX?
bool cHardwareBCR::Inst_SetOpinion(cAvidaContext& ctx)
{
  assert(m_organism != 0);
  m_organism->GetOrgInterface().SetOpinion(GetRegister(FindModifiedRegister(rBX)), m_organism);
  return true;
}

/* Gets the organism's current opinion, placing the opinion in register ?BX?
   and the age of the opinion in register !?BX?
 */
bool cHardwareBCR::Inst_GetOpinion(cAvidaContext& ctx)
{
  assert(m_organism != 0);
  if (m_organism->GetOrgInterface().HasOpinion(m_organism)) {
    const int opinion_reg = FindModifiedRegister(rBX);
    const int age_reg = FindNextRegister(opinion_reg);

    setInternalValue(opinion_reg, m_organism->GetOpinion().first, true);
    setInternalValue(age_reg, m_world->GetStats().GetUpdate() - m_organism->GetOpinion().second, true);
  }
  return true;
}

//! An organism joins a group by setting it opinion to the group id. 
bool cHardwareBCR::Inst_JoinGroup(cAvidaContext&)
{
  int opinion = m_world->GetConfig().DEFAULT_GROUP.Get();
  // Check if the org is currently part of a group
  assert(m_organism != 0);
	
  int prop_group_id = GetRegister(FindModifiedRegister(rBX));
  
  // check if this is a valid group
  if (m_world->GetConfig().USE_FORM_GROUPS.Get() == 2 &&
      !(m_world->GetEnvironment().IsGroupID(prop_group_id))) {
    return false; 
  }
  // injected orgs might not have an opinion
  if (m_organism->GetOrgInterface().HasOpinion(m_organism)) {
    opinion = m_organism->GetOpinion().first;
    
    //return false if org setting opinion to current one (avoid paying costs for not switching)
    if (opinion == prop_group_id) return false;
    
    // A random chance for failure to join group based on config, if failed return true for resource cost.
    if (m_world->GetConfig().JOIN_GROUP_FAILURE.Get() > 0) {
      int percent_failure = m_world->GetConfig().JOIN_GROUP_FAILURE.Get();
      double prob_failure = (double) percent_failure / 100.0;
      double rand = m_world->GetRandom().GetDouble();
      if (rand <= prob_failure) return true;
    }
    
    // If tolerances are on the org must pass immigration chance 
    if (m_world->GetConfig().TOLERANCE_WINDOW.Get() > 0) {
      m_organism->GetOrgInterface().AttemptImmigrateGroup(prop_group_id, m_organism);
      return true;
    }
    else {
      // otherwise, subtract org from current group
      m_organism->LeaveGroup(opinion);
    }
  }
  
  // Set the opinion
  m_organism->GetOrgInterface().SetOpinion(prop_group_id, m_organism);
  
  // Add org to group count
  if (m_organism->GetOrgInterface().HasOpinion(m_organism)) {
    opinion = m_organism->GetOpinion().first;	
    m_organism->JoinGroup(opinion);
  }
  
  return true;
}

bool cHardwareBCR::Inst_GetGroupID(cAvidaContext&)
{
  assert(m_organism != 0);
  if (m_organism->HasOpinion()) {
    const int group_reg = FindModifiedRegister(rBX);
    
    setInternalValue(group_reg, m_organism->GetOpinion().first, false);
  }
  return true;
}

bool cHardwareBCR::Inst_GetFacedOrgID(cAvidaContext&)
//Get ID of organism faced by this one, if there is an organism in front.
{
  if (m_use_avatar && m_use_avatar != 2) return false;
  cOrganism* neighbor = NULL;
  if (!m_use_avatar && !m_organism->IsNeighborCellOccupied()) return false;
  else if (m_use_avatar == 2 && !m_organism->GetOrgInterface().FacedHasAV()) return false;
  
  if (!m_use_avatar) neighbor = m_organism->GetOrgInterface().GetNeighbor();
  else if (m_use_avatar == 2) neighbor = m_organism->GetOrgInterface().GetRandFacedAV();
  if (neighbor->IsDead())  return false;  
  
  const int out_reg = FindModifiedRegister(rBX);
  setInternalValue(out_reg, neighbor->GetID(), true);
  return true;
}

//Teach offspring learned targeting/foraging behavior
bool cHardwareBCR::Inst_TeachOffspring(cAvidaContext&)
{
  assert(m_organism != 0);
  m_organism->Teach(true);
  return true;
}

bool cHardwareBCR::Inst_LearnParent(cAvidaContext& ctx)
{
  assert(m_organism != 0);
  bool halt = false;
  if (m_organism->HadParentTeacher()) {
    int old_target = m_organism->GetForageTarget();
    int prop_target = -1;
    prop_target = m_organism->GetParentFT();

    halt = (prop_target <= -2 && m_world->GetConfig().PRED_PREY_SWITCH.Get() < 0);
    if (!halt) {
      if (m_use_avatar && m_organism->GetOrgInterface().GetAVCellID() != -1 && 
          (((prop_target == -2 || prop_target == -3) && old_target > -2) || (prop_target > -2 && (old_target == -2 || prop_target == -3)))) {
        m_organism->GetOrgInterface().SwitchPredPrey();
        m_organism->CopyParentFT(ctx);
      }
      else m_organism->CopyParentFT(ctx);
    }
  }
  return !halt;
}

bool cHardwareBCR::Inst_ScrambleReg(cAvidaContext& ctx)
{
  for (int i = 0; i < NUM_REGISTERS; i++) {
    setInternalValue(rAX + i, (int) (ctx.GetRandom().GetDouble()), true);
  }
  return true;
}

bool cHardwareBCR::DoActualCollect(cAvidaContext& ctx, int bin_used, bool unit)
{
  // Set up res_change and max total
  Apto::Array<double> res_count;
  if (!m_use_avatar) res_count = m_organism->GetOrgInterface().GetResources(ctx);
  else if (m_use_avatar) res_count = m_organism->GetOrgInterface().GetAVResources(ctx); 
  Apto::Array<double> res_change(res_count.GetSize());
  res_change.SetAll(0.0);
  double total = m_organism->GetRBinsTotal();
  double max = m_world->GetConfig().MAX_TOTAL_STORED.Get();
  bool has_max = max > 0 ? true : false;
  double res_consumed = 0.0;
  
  // Collect a unit or some ABSORB_RESOURCE_FRACTION
  const cResourceLib& resource_lib = m_world->GetEnvironment().GetResourceLib();
  if (unit) {
    double threshold = resource_lib.GetResource(bin_used)->GetThreshold();
    if (res_count[bin_used] >= threshold) {
      res_consumed = threshold;
    }
    else {
      return false;
    }  
  }
  else {
    res_consumed = res_count[bin_used] * m_world->GetConfig().ABSORB_RESOURCE_FRACTION.Get();
  }
  
  if (has_max && res_consumed + total >= max) {
    res_consumed = max - total;
    res_change[bin_used] = -1 * res_consumed;
  }
  else res_change[bin_used] = -1 * res_consumed;
  
  if (res_consumed > 0) { 
    m_organism->AddToRBin(bin_used, res_consumed); 
    if (!m_use_avatar) m_organism->GetOrgInterface().UpdateResources(ctx, res_change);
    else if (m_use_avatar) m_organism->GetOrgInterface().UpdateAVResources(ctx, res_change);
    return true;
  }
  return false;
}