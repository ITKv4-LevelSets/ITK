/*=========================================================================

  Program:   Insight Segmentation & Registration Toolkit
  Module:    itkWatershedSegmentTreeGenerator.txx
  Language:  C++
  Date:      $Date$
  Version:   $Revision$

  Copyright (c) 2002 Insight Consortium. All rights reserved.
  See ITKCopyright.txt or http://www.itk.org/HTML/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even 
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR 
     PURPOSE.  See the above copyright notices for more information.

=========================================================================*/
#ifndef __itkWatershedSegmentTreeGenerator_txx
#define __itkWatershedSegmentTreeGenerator_txx

#include <stack>
#include "itkWatershedOneWayEquivalencyTable.h"

namespace itk
{
namespace watershed
{

template <class TScalarType>
SegmentTreeGenerator<TScalarType>
::SegmentTreeGenerator() : m_Merge(false), m_FloodLevel(0.0),
  m_ConsumeInput(false)
{
  typename SegmentTreeType::Pointer st
    = static_cast<SegmentTreeType*>(this->MakeOutput(0).GetPointer());
  this->SetNumberOfRequiredOutputs(1);
  this->ProcessObject::SetNthOutput(0, st.GetPointer());
  m_MergedSegmentsTable = OneWayEquivalencyTableType::New();
}

template <class TScalarType>
typename SegmentTreeGenerator<TScalarType>::DataObjectPointer
SegmentTreeGenerator<TScalarType>
::MakeOutput(unsigned int idx)
{
  return static_cast<DataObject*>(SegmentTreeType::New().GetPointer());
}
  
template <class TScalarType>
void SegmentTreeGenerator<TScalarType>
::SetFloodLevel(double val)
{
  // Clamp level between 0.0 and 1.0
  if (val > 1.0)         { m_FloodLevel = 1.0; }
  else if (val < 0.0)    { m_FloodLevel = 0.0; }
  else                   { m_FloodLevel = val; }

  // Has this flood level increased over a level that has already been
  // calculated? If not, then we don't set the modified flag.  A decrease in
  // the flood level does not require re-execution of the filter.
  if ( m_HighestCalculatedFloodLevel < m_FloodLevel )
    {
      this->SetHighestCalculatedFloodLevel( m_FloodLevel );
      this->Modified(); // redundant
    }
}
  
template <class TScalarType>
void SegmentTreeGenerator<TScalarType>
::GenerateData()
{
  m_MergedSegmentsTable->Clear();

  typename SegmentTableType::Pointer input = this->GetInputSegmentTable();
  typename SegmentTreeType::Pointer mergeList   = SegmentTreeType::New();
  typename SegmentTableType::Pointer seg = SegmentTableType::New();
  if (m_ConsumeInput == true) // do not copy input
    {
      input->Modified();
      input->SortEdgeLists();
      if (m_Merge == true)   {      this->MergeEquivalencies();    }
      this->CompileMergeList(input, mergeList);
      this->ExtractMergeHierarchy(input, mergeList);
    }
  else
    {
      seg->Copy(*input); // copy the input
      seg->SortEdgeLists();
      if (m_Merge == true)   {      this->MergeEquivalencies();    }
      this->CompileMergeList(seg, mergeList);
      this->ExtractMergeHierarchy(seg, mergeList);
    }
  this->UpdateProgress(1.0);
}

template <class TScalarType>
void SegmentTreeGenerator<TScalarType>
::MergeEquivalencies()
{
  typename SegmentTableType::Pointer segTable = this->GetInputSegmentTable();
  typename EquivalencyTableType::Pointer eqTable  =
    this->GetInputEquivalencyTable();
  typename EquivalencyTableType::Iterator it;
  ScalarType threshold = m_FloodLevel * segTable->GetMaximumDepth();

  eqTable->Flatten();
  unsigned long counter =0;

  segTable->PruneEdgeLists(threshold);

  for (it = eqTable->Begin(); it != eqTable->End(); ++it)
    {
      MergeSegments(segTable, m_MergedSegmentsTable,
            (*it).first, (*it).second);  // Merge first INTO second.
                                         // deletes first
      if ( (counter % 10000) == 0 )
        {
          segTable->PruneEdgeLists(threshold);
          m_MergedSegmentsTable->Flatten();
          counter = 0;
        }

      counter++;
    }
}

template <class TScalarType>
void SegmentTreeGenerator<TScalarType>
::CompileMergeList(SegmentTableTypePointer segments,
                   SegmentTreeTypePointer mergeList)
{
  typename SegmentTreeType::merge_t tempMerge;

  // Region A will flood Region B (B will merge with A) at a flood level L
  // when all of the following conditions are true:
  // 1) Depth of B < L
  // 2) A is across the lowest edge of B
  typename SegmentTableType::Iterator segment_ptr;
  unsigned long labelFROM;
  unsigned long labelTO;
  ScalarType threshold = m_FloodLevel * segments->GetMaximumDepth();
  m_MergedSegmentsTable->Flatten();

  segments->PruneEdgeLists(threshold);
  
  for (segment_ptr = segments->Begin(); segment_ptr != segments->End();
       ++segment_ptr)
    {
      labelFROM = (*segment_ptr).first;
      
      // Must take into account any equivalencies that have already been
      // recorded.
      labelTO
        = m_MergedSegmentsTable->RecursiveLookup((*segment_ptr).second.edge_list.front().label);
      while (labelTO == labelFROM) // Pop off any bogus merges with ourself
        {                          // that may have been left in this list.
          (*segment_ptr).second.edge_list.pop_front();
          labelTO
            = m_MergedSegmentsTable->RecursiveLookup((*segment_ptr).second.edge_list.front().label);
        }

      // Add this merge to our list if its saliency is below
      // the threshold.
      tempMerge.from     = labelFROM;
      tempMerge.to       = labelTO;
      tempMerge.saliency = (*segment_ptr).second.edge_list.front().height
        - (*segment_ptr).second.min;
      if (tempMerge.saliency < threshold)
        {
          mergeList->PushBack(tempMerge);
        }
    }
  
  // Heapsort the list
  typename SegmentTreeType::merge_comp comp;
  std::make_heap(mergeList->Begin(), mergeList->End(), comp);
}

template <class TScalarType>
void SegmentTreeGenerator<TScalarType>
::ExtractMergeHierarchy(SegmentTableTypePointer segments,
                        SegmentTreeTypePointer heap)
{
  typename SegmentTreeType::Pointer list = this->GetOutputSegmentTree();

  // Merges segments up to a specified floodlevel according to the information
  // in the heap of merges.  As two segments are merged, calculates a new
  // possible merges and pushes it onto the heap.
  ScalarType threshold = m_FloodLevel * segments->GetMaximumDepth();

  unsigned  counter;
  typename SegmentTreeType::merge_comp comp;
  typename SegmentTableType::DataType  *toSeg;
  typename SegmentTreeType::ValueType  tempMerge;
  unsigned long  toSegLabel, fromSegLabel;

  if (heap->Empty()) return;
  double initHeapSize = static_cast<double>(heap->Size());

  counter = 0;
  typename SegmentTreeType::ValueType topMerge = heap->Front();

  while( (! heap->Empty()) && (topMerge.saliency <= threshold) )
    {
      counter++;             // Every so often we should eliminate
      if ((counter == 10000)) // all the recursion in our records
        {                    // of which segments have merged.
          counter = 0;
          segments->PruneEdgeLists(threshold); // also we want to
          // keep the edge list size under control
        }
      if ((counter % 10000) == 0)
        {
          m_MergedSegmentsTable->Flatten();
        }

      if ((counter % 1000) == 0)
        {
          this->UpdateProgress(1.0 - ((static_cast<double>(heap->Size())) /
                                      initHeapSize));
        }
      
      std::pop_heap(heap->Begin(), heap->End(), comp);
      heap->PopBack();  // Popping the heap moves the top element to the end
                        // of the container structure, so we delete that here.

      // Recursively find the segments we are about to merge
      // (the labels identified here may have merged already)
      fromSegLabel = m_MergedSegmentsTable->RecursiveLookup(topMerge.from);
      toSegLabel   = m_MergedSegmentsTable->RecursiveLookup(topMerge.to);

      // If the two segments do not resolve to the same segment and the
      // "TO" segment has never been merged, then then merge them.
      // Otherwise, ignore this particular entry.
      if ( fromSegLabel == topMerge.from && fromSegLabel != toSegLabel )
        {
          toSeg   = segments->Lookup( toSegLabel );

          topMerge.from = fromSegLabel;
          topMerge.to   = toSegLabel;
          list->PushBack(topMerge);  // Record this merge for posterity

          // Merge the segments
          Self::MergeSegments(segments, m_MergedSegmentsTable,
                              fromSegLabel, toSegLabel);

          // Now check for new possible merges in A.
          // All we have to do is look at the front of the
          // ordered list.

          // Recursively look up the label to which we might
          // be merging to.
          if (! toSeg->edge_list.empty() )
            {
              tempMerge.from = toSegLabel;  // The new, composite segment
              tempMerge.to   = m_MergedSegmentsTable->RecursiveLookup(
                                        toSeg->edge_list.front().label );
              while (tempMerge.to == tempMerge.from)
                { // We don't want to merge to ourself.
                  toSeg->edge_list.pop_front();
                  tempMerge.to = m_MergedSegmentsTable->RecursiveLookup(
                                        toSeg->edge_list.front().label );
                }
              tempMerge.saliency =
                (toSeg->edge_list.front().height) - toSeg->min;
              
              heap->PushBack(tempMerge);
              std::push_heap(heap->Begin(), heap->End(), comp);
            }
        }
      if( ! heap->Empty() )
        {
          topMerge = heap->Front();
        }
    }
}


template <class TScalarType>
void SegmentTreeGenerator<TScalarType>
::PruneMergeSegments(SegmentTableTypePointer segments,
                     OneWayEquivalencyTableTypePointer eqT,
                     const unsigned long FROM, const unsigned long TO,
                     ScalarType maxSaliency)
{
  typename SegmentTableType::edge_list_t::iterator edgeTOi, edgeFROMi,
    edgeTEMPi; 
  itk::hash_map<unsigned long, bool, itk::hash<unsigned long> >
    seen_table;
  unsigned long labelTO, labelFROM;

  // Lookup both entries.
  typename SegmentTableType::segment_t *from_seg = segments->Lookup(FROM);
  typename SegmentTableType::segment_t *to_seg   = segments->Lookup(TO);

  if (from_seg == 0 || to_seg == 0)
    {
      ExceptionObject e(__FILE__, __LINE__);
      std::ostrstream msg;
      msg << "SegmentTreeGenerator::PruneMergeSegments()" << std::ends;
      e.SetLocation(msg.str());
      e.SetDescription("An unexpected and fatal error has occurred.");
      throw e;
    }

  // Compare the minimum values.
  if ((from_seg->min) < (to_seg->min))  (to_seg->min) = (from_seg->min);

  // Add all the "FROM" edges to the "TO" edges.
  
  // We MUST eliminate redundant edges to prevent the
  // segment table from growing to unmanageable sizes.
  // Since the lists are sorted, the lowest height
  // value will be added first, so any redundant edges
  // can simply be eliminated.

  // References in adjacent edge tables are not replaced,
  // but rather will be resolved later through the one-way
  // equivalency table.
  
  edgeTOi   = to_seg->edge_list.begin();
  edgeFROMi = from_seg->edge_list.begin();
  while ( edgeTOi != to_seg->edge_list.end() &&
          edgeFROMi != from_seg->edge_list.end() )
    {
      // Recursively resolve the labels we are seeing
      labelTO   = eqT->RecursiveLookup(edgeTOi->label);
      labelFROM = eqT->RecursiveLookup(edgeFROMi->label);

      // Ignore any labels already in this list and
      // any pointers back to ourself.
      // This step is necessary to keep the edge lists from
      // growing exponentially in size as segments merge.
      if ( seen_table.find(labelTO) != seen_table.end() ||
           labelTO == FROM )
        {
          edgeTEMPi = edgeTOi;
          edgeTEMPi++;
          to_seg->edge_list.erase(edgeTOi);
          edgeTOi = edgeTEMPi;
          continue;
        }
      if ( seen_table.find(labelFROM) != seen_table.end() ||
           labelFROM == TO)
        {
          edgeFROMi++;
          continue;
        }

      // Fix any changed labels since we have the chance
      if ( labelTO   != edgeTOi->label  ) edgeTOi->label = labelTO;
      if ( labelFROM != edgeFROMi->label) edgeFROMi->label = labelFROM;

      // Which edge is next in the list?
      if ( edgeFROMi->height < edgeTOi->height )
        {
          to_seg->edge_list.insert(edgeTOi, *edgeFROMi);
          seen_table.insert(std::pair<unsigned long, bool>(labelFROM, true));
          edgeFROMi++;
        }
      else
        {
          seen_table.insert(std::pair<unsigned long, bool>(labelTO, true));
          edgeTOi++;
        }
    }
  
  // Process tail of the FROM list.
  while ( edgeFROMi != from_seg->edge_list.end() )
    {
      labelFROM = eqT->RecursiveLookup(edgeFROMi->label);
      if ( seen_table.find(labelFROM) != seen_table.end() ||
           labelFROM == TO)
        {
          edgeFROMi++;
        }
      else
        {
          if ( labelFROM != edgeFROMi->label) edgeFROMi->label = labelFROM;
          to_seg->edge_list.push_back(*edgeFROMi);
          seen_table.insert(std::pair<unsigned long, bool>(labelFROM, true));
          edgeFROMi++;
        }
    }

  // Process tail of the TO list.
  while ( edgeTOi != to_seg->edge_list.end() )
    {
      labelTO   = eqT->RecursiveLookup(edgeTOi->label);
      if ( seen_table.find(labelTO) != seen_table.end() ||
           labelTO == FROM)
        {
          edgeTEMPi = edgeTOi;
          edgeTEMPi++;
          to_seg->edge_list.erase(edgeTOi);
          edgeTOi = edgeTEMPi;
        }
      else
        {
          if ( labelTO   != edgeTOi->label  ) edgeTOi->label = labelTO;
          seen_table.insert(std::pair<unsigned long, bool>(labelTO, true));
          edgeTOi++;
        }
    }

  // Now destroy the segment that was merged.
  segments->Erase(FROM);

  // Record this equivalency
  eqT->Add(FROM, TO);
}

template <class TScalarType>
void SegmentTreeGenerator<TScalarType>
::MergeSegments(SegmentTableTypePointer segments,
                OneWayEquivalencyTableTypePointer eqT,
                const unsigned long FROM, const unsigned long TO)
{
  typename SegmentTableType::edge_list_t::iterator edgeTOi, edgeFROMi,
    edgeTEMPi;
  itk::hash_map<unsigned long, bool, itk::hash<unsigned long> >
    seen_table;
  unsigned long labelTO, labelFROM;

  // Lookup both entries.
  typename SegmentTableType::segment_t *from_seg = segments->Lookup(FROM);
  typename SegmentTableType::segment_t *to_seg   = segments->Lookup(TO);

  if (from_seg == 0 || to_seg == 0)
    {
      ExceptionObject e(__FILE__, __LINE__);
      std::ostrstream msg;
      msg <<"SegmentTreeGenerator::MergeSegments()" << std::ends;
      e.SetLocation(msg.str());
      e.SetDescription("An unexpected and fatal error has occurred.");
      throw e;
    }

  // Compare the minimum values.
  if ((from_seg->min) < (to_seg->min))  (to_seg->min) = (from_seg->min);

  // Add all the "FROM" edges to the "TO" edges.
  
  // We MUST eliminate redundant edges to prevent the
  // segment table from growing to unmanageable sizes.
  // Since the lists are sorted, the lowest height
  // value will be added first, so any redundant edges
  // can simply be eliminated.

  // References in adjacent edge tables are not replaced,
  // but rather will be resolved later through the one-way
  // equivalency table.
  
  edgeTOi   = to_seg->edge_list.begin();
  edgeFROMi = from_seg->edge_list.begin();
  while ( edgeTOi != to_seg->edge_list.end() &&
          edgeFROMi != from_seg->edge_list.end() )
    {
      // Recursively resolve the labels we are seeing
      labelTO   = eqT->RecursiveLookup(edgeTOi->label);
      labelFROM = eqT->RecursiveLookup(edgeFROMi->label);

      // Ignore any labels already in this list and
      // any pointers back to ourself.
      // This step is necessary to keep the edge lists from
      // growing exponentially in size as segments merge.
      if ( seen_table.find(labelTO) != seen_table.end() ||
           labelTO == FROM )
        {
          edgeTEMPi = edgeTOi;
          edgeTEMPi++;
          to_seg->edge_list.erase(edgeTOi);
          edgeTOi = edgeTEMPi;
          continue;
        }
      if ( seen_table.find(labelFROM) != seen_table.end() ||
           labelFROM == TO)
        {
          edgeFROMi++;
          continue;
        }

      // Fix any changed labels since we have the chance
      if ( labelTO   != edgeTOi->label  ) edgeTOi->label = labelTO;
      if ( labelFROM != edgeFROMi->label) edgeFROMi->label = labelFROM;

      // Which edge is next in the list?
      if ( edgeFROMi->height < edgeTOi->height )
        {
          to_seg->edge_list.insert(edgeTOi, *edgeFROMi);
          seen_table.insert(std::pair<unsigned long, bool>(labelFROM, true));
          edgeFROMi++;
        }
      else
        {
          seen_table.insert(std::pair<unsigned long, bool>(labelTO, true));
          edgeTOi++;
        }
    }
  
  // Process tail of the FROM list.
  while ( edgeFROMi != from_seg->edge_list.end() )
    {
      labelFROM = eqT->RecursiveLookup(edgeFROMi->label);
      if ( seen_table.find(labelFROM) != seen_table.end() ||
           labelFROM == TO)
        {
          edgeFROMi++;
        }
      else
        {
          if ( labelFROM != edgeFROMi->label) edgeFROMi->label = labelFROM;
          to_seg->edge_list.push_back(*edgeFROMi);
          seen_table.insert(std::pair<unsigned long, bool>(labelFROM, true));
          edgeFROMi++;
        }
    }

  // Process tail of the TO list.
  while ( edgeTOi != to_seg->edge_list.end() )
    {
      labelTO   = eqT->RecursiveLookup(edgeTOi->label);
      if ( seen_table.find(labelTO) != seen_table.end() ||
           labelTO == FROM)
        {
          edgeTEMPi = edgeTOi;
          edgeTEMPi++;
          to_seg->edge_list.erase(edgeTOi);
          edgeTOi = edgeTEMPi;
        }
      else
        {
          if ( labelTO   != edgeTOi->label  ) edgeTOi->label = labelTO;
          seen_table.insert(std::pair<unsigned long, bool>(labelTO, true));
          edgeTOi++;
        }
    }

  // Now destroy the segment that was merged.
  segments->Erase(FROM);

  // Record this equivalency
  eqT->Add(FROM, TO);  
}

template <class TScalarType>  
void SegmentTreeGenerator<TScalarType>
::GenerateOutputRequestedRegion(DataObject *output)
{
}

template <class TScalarType>  
void SegmentTreeGenerator<TScalarType>
::GenerateInputRequestedRegion()
{
  Superclass::GenerateInputRequestedRegion();

   // get pointers to the input and output
  //  typename SegmentTableType::Pointer  inputPtr  = this->GetInputSegmentTable();
  //  typename SegmentTreeType::Pointer outputPtr = this->GetOutputSegmentTree();
  
  //  if ( !inputPtr || !outputPtr )
  //    {
  //    return;
  //    }

}

template<class TScalarType>
void 
SegmentTreeGenerator<TScalarType>
::PrintSelf(std::ostream& os, Indent indent) const
{
  Superclass::PrintSelf(os,indent);
  os << indent << "FloodLevel: " << m_FloodLevel << std::endl;
  os << indent << "Merge: " << m_Merge << std::endl;
  os << indent << "ConsumeInput: " << m_ConsumeInput << std::endl;
  os << indent << "HighestCalculatedFloodLevel: " <<
    m_HighestCalculatedFloodLevel << std::endl;
}  
 
}// end namespace watershed
}// end namespace itk

#endif
