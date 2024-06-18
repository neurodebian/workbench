/*LICENSE_START*/
/*
 *  Copyright (C) 2014  Washington University School of Medicine
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License along
 *  with this program; if not, write to the Free Software Foundation, Inc.,
 *  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
 */
/*LICENSE_END*/

#include "AlgorithmCiftiReplaceStructure.h"
//for computing the cropped volume space
#include "AlgorithmCiftiSeparate.h"
#include "AlgorithmException.h"
#include "CaretLogger.h"
#include "CaretPointer.h"
#include "CiftiFile.h"
#include "GiftiLabelTable.h"
#include "LabelFile.h"
#include "MetricFile.h"
#include "Vector3D.h"
#include "VolumeFile.h"

#include <cmath>
#include <map>
#include <set>

using namespace caret;
using namespace std;

AString AlgorithmCiftiReplaceStructure::getCommandSwitch()
{
    return "-cifti-replace-structure";
}

AString AlgorithmCiftiReplaceStructure::getShortDescription()
{
    return "REPLACE DATA IN A STRUCTURE IN A CIFTI FILE";
}

OperationParameters* AlgorithmCiftiReplaceStructure::getParameters()
{
    OperationParameters* ret = new OperationParameters();
    ret->addStringParameter(1, "cifti", "the cifti to modify");//see useParameters for why this is a string
    
    ret->addStringParameter(2, "direction", "which dimension to interpret as a single map, ROW or COLUMN");
    
    ParameterComponent* labelOpt = ret->createRepeatableParameter(3, "-label", "replace the data in a surface label component");
    labelOpt->addStringParameter(1, "structure", "the structure to replace the data of");
    labelOpt->addLabelParameter(2, "label", "the input label file");

    ParameterComponent* metricOpt = ret->createRepeatableParameter(4, "-metric", "replace the data in a surface component");
    metricOpt->addStringParameter(1, "structure", "the structure to replace the data of");
    metricOpt->addMetricParameter(2, "metric", "the input metric");

    ParameterComponent* volumeOpt = ret->createRepeatableParameter(5, "-volume", "replace the data in a volume component");
    volumeOpt->addStringParameter(1, "structure", "the structure to replace the data of");
    volumeOpt->addVolumeParameter(2, "volume", "the input volume");
    volumeOpt->createOptionalParameter(3, "-from-cropped", "the input is cropped to the size of the component");

    OptionalParameter* volumeAllOpt = ret->createOptionalParameter(6, "-volume-all", "replace the data in all volume components");
    volumeAllOpt->addVolumeParameter(1, "volume", "the input volume");
    volumeAllOpt->createOptionalParameter(2, "-from-cropped", "the input is cropped to the size of the data");
    
    ret->createOptionalParameter(7, "-discard-unused-labels", "when operating on a dlabel file, drop any unused label keys from the label table");
    
    OptionalParameter* collideOpt = ret->createOptionalParameter(8, "-label-collision", "how to handle conflicts between label keys");
    collideOpt->addStringParameter(1, "action", "'ERROR', 'LEFT_SURFACE_FIRST', or 'LEGACY', default 'ERROR', use 'LEGACY' to match v1.4.2 and earlier");

    AString helpText = AString("This is a fairly low-level command, you probably want to use -cifti-create-dense-from-template instead.\n\n") +
        "You must specify at least one of -metric, -label, -volume, or -volume-all for this command to do anything.  " +
        "Input volumes must line up with the output of -cifti-separate.  " +
        "For dtseries/dscalar, use COLUMN, and if your dconn matrix will be fully symmetric, COLUMN is more efficient.  " +
        "The -volume-all option must not be specified when using a -volume option.  " +
        "A -metric option must not be specified when using a -label option, and is not recommended on a label-type cifti file.  " +
        "For each <structure> argument, use one of the following strings:\n";
    vector<StructureEnum::Enum> myStructureEnums;
    StructureEnum::getAllEnums(myStructureEnums);
    for (int i = 0; i < (int)myStructureEnums.size(); ++i)
    {
        helpText += "\n" + StructureEnum::toName(myStructureEnums[i]);
    }
    ret->setHelpText(helpText);
    return ret;
}

namespace
{
    enum DataSourceType
    {
        LABEL,
        METRIC,
        VOLUME,
        VOLUME_ALL
    };
    
    enum LabelConflictLogic
    {
        ERROR,
        LEFT_SURFACE_FIRST,
        LEGACY
    };

    struct DataSourceInfo
    {
        DataSourceType type;
        int64_t index;//for repeatable options

        DataSourceInfo(const DataSourceType myType, const int64_t myIndex)
        {
            type = myType;
            index = myIndex;
        };
        DataSourceInfo()
        {
            index = -1;//make it obviously invalid by default
        };
    };
}

void AlgorithmCiftiReplaceStructure::useParameters(OperationParameters* myParams, ProgressObject* /*myProgObj*/)
{
    AString ciftiName = myParams->getString(1);
    CiftiFile myCifti(ciftiName);//CiftiFile currently doesn't expose the filename it is using, and we need to call writeFile manually (since the parsing framework doesn't expect an input to be modified)
    AString dirName = myParams->getString(2);//as for why it needs to modify an input, that is the only way it is useful for things like -cifti-smoothing
    int myDir;
    if (dirName == "ROW")
    {
        myDir = CiftiXML::ALONG_ROW;
    } else if (dirName == "COLUMN") {
        myDir = CiftiXML::ALONG_COLUMN;
    } else {
        throw AlgorithmException("incorrect string for direction, use ROW or COLUMN");
    }
    bool discardUnusedLabels = myParams->getOptionalParameter(7)->m_present;
    LabelConflictLogic conflictLogic = ERROR;
    OptionalParameter* collideOpt = myParams->getOptionalParameter(8);
    if (collideOpt->m_present)
    {
        AString collideStr = collideOpt->getString(1);
        if (collideStr == "ERROR")
        {
            conflictLogic = ERROR;
        } else if (collideStr == "LEFT_SURFACE_FIRST") {
            conflictLogic = LEFT_SURFACE_FIRST;
        } else if (collideStr == "LEGACY") {
            conflictLogic = LEGACY;
        } else {
            throw AlgorithmException("incorrect string for label collision option");
        }
    }
    const vector<ParameterComponent*>& labelInst = myParams->getRepeatableParameterInstances(3);
    const vector<ParameterComponent*>& metricInst = myParams->getRepeatableParameterInstances(4);
    const vector<ParameterComponent*>& volumeInst = myParams->getRepeatableParameterInstances(5);
    OptionalParameter* volumeAllOpt = myParams->getOptionalParameter(6);
    if (volumeAllOpt->m_present && volumeInst.size() != 0)
    {
        throw AlgorithmException("using -volume-all and -volume at the same time in -cifti-replace-structure is strongly discouraged");
    }
    map<StructureEnum::Enum, DataSourceInfo> surfSources, volSources;
    vector<DataSourceInfo> legacyOrder;//label versus metric was NOT checked for metric replace on dlabel...but, we should probably still error when the arguments are very wrong
    int labelMode = 0;//0 unknown, 1 label, 2 not label
    for (int i = 0; i < (int)labelInst.size(); ++i)
    {
        labelMode = 1;
        AString structName = labelInst[i]->getString(1);
        bool ok = false;
        StructureEnum::Enum myStruct = StructureEnum::fromName(structName, &ok);
        if (!ok)
        {
            throw AlgorithmException("unrecognized structure name");
        }
        if (surfSources.find(myStruct) != surfSources.end())
        {
            throw AlgorithmException("-label specified more than once with structure " + structName);//this wasn't previously checked either...
        }
        DataSourceInfo thisInfo(LABEL, i);
        surfSources[myStruct] = thisInfo;
        legacyOrder.push_back(thisInfo);
    }
    for (int i = 0; i < (int)metricInst.size(); ++i)
    {
        if (labelMode == 1)
        {
            throw AlgorithmException("-metric and -label options cannot be used together");
        }
        labelMode = 2;
        AString structName = metricInst[i]->getString(1);
        bool ok = false;
        StructureEnum::Enum myStruct = StructureEnum::fromName(structName, &ok);
        if (!ok)
        {
            throw AlgorithmException("unrecognized structure name");
        }
        if (surfSources.find(myStruct) != surfSources.end())
        {
            throw AlgorithmException("-metric specified more than once with structure " + structName);//this wasn't previously checked either...
        }
        DataSourceInfo thisInfo(METRIC, i);
        surfSources[myStruct] = thisInfo;
        legacyOrder.push_back(thisInfo);
    }
    for (int i = 0; i < (int)volumeInst.size(); ++i)
    {
        AString structName = volumeInst[i]->getString(1);
        bool ok = false;
        StructureEnum::Enum myStruct = StructureEnum::fromName(structName, &ok);
        if (!ok)
        {
            throw AlgorithmException("unrecognized structure name");
        }
        if (volSources.find(myStruct) != volSources.end())
        {
            throw AlgorithmException("-volume specified more than once with structure " + structName);//this wasn't previously checked either...
        }
        DataSourceInfo thisInfo(VOLUME, i);
        volSources[myStruct] = thisInfo;
        legacyOrder.push_back(thisInfo);
        bool thisLabel = (volumeInst[i]->getVolume(2)->getType() == SubvolumeAttributes::LABEL);
        int thisLabelMode = (thisLabel ? 1 : 2);
        if (labelMode == 0)
        {
            labelMode = thisLabelMode;
        } else {
            if (labelMode != thisLabelMode)
            {
                if (thisLabel)
                {
                    throw AlgorithmException("label volume encountered when other inputs are non-label type: " + volumeInst[i]->getVolume(2)->getFileName());
                } else {
                    throw AlgorithmException("non-label volume encountered when other inputs are label type: " + volumeInst[i]->getVolume(2)->getFileName());
                }
            }
        }
    }
    if (volumeAllOpt->m_present)
    {
        legacyOrder.push_back(DataSourceInfo(VOLUME_ALL, -1));
        bool thisLabel = (volumeAllOpt->getVolume(1)->getType() == SubvolumeAttributes::LABEL);
        int thisLabelMode = (thisLabel ? 1 : 2);
        if (labelMode == 0)
        {
            labelMode = thisLabelMode;
        } else {
            if (labelMode != thisLabelMode)
            {
                if (thisLabel)
                {
                    throw AlgorithmException("label volume encountered when other inputs are non-label type: " + volumeAllOpt->getVolume(1)->getFileName());
                } else {
                    throw AlgorithmException("non-label volume encountered when other inputs are label type: " + volumeAllOpt->getVolume(1)->getFileName());
                }
            }
        }
    }
    bool errorOnLabelConflict = (legacyOrder.size() > 1) && (conflictLogic == ERROR);//don't error when all priority options would give the same result
    if (errorOnLabelConflict && labelMode == 1)//check for replacing everything when the new label tables don't collide, as that behavior also doesn't depend on priority
    {
        const CiftiXML& myXML = myCifti.getCiftiXML();
        if (myXML.getMappingType(myDir) != CiftiMappingType::BRAIN_MODELS) throw AlgorithmException("specified dimension must contain brain models");
        const CiftiBrainModelsMap& myMap = myXML.getBrainModelsMap(myDir);
        bool allReplace = true;
        for (auto surfStruct : myMap.getSurfaceStructureList())
        {
            if (surfSources.find(surfStruct) == surfSources.end())
            {
                allReplace = false;
                break;
            }
        }
        for (auto volStruct : myMap.getVolumeStructureList())
        {
            if (volSources.find(volStruct) == volSources.end())
            {
                allReplace = false;
                break;
            }
        }
        if (allReplace)
        {
            errorOnLabelConflict = false;//replace everything means we can ignore conflicts with the labels in the original cifti
            vector<GiftiLabelTable> testTables;//but need to manually check the tables of the inputs
            for (auto source : legacyOrder)//order doesn't really matter if we are just generating an error
            {
                CaretMappableDataFile* thisFile = NULL;
                switch (source.type)
                {
                    case LABEL:
                        thisFile = labelInst[source.index]->getLabel(2);
                        break;
                    case VOLUME:
                        thisFile = volumeInst[source.index]->getVolume(2);
                        break;
                    case VOLUME_ALL:
                        thisFile = volumeAllOpt->getVolume(1);
                        break;
                    default:
                        CaretAssert(false);//it is an error to mix -metric and -label now, so this can't happen
                        throw AlgorithmException("internal error, tell the developers what you tried to do");
                }
                if (testTables.empty())
                {
                    testTables.resize(thisFile->getNumberOfMaps());
                } else {
                    if (thisFile->getNumberOfMaps() != int(testTables.size()))
                    {
                        throw AlgorithmException("input files have different numbers of maps");
                    }
                }
                for (int i = 0; i < int(testTables.size()); ++i)
                {
                    testTables[i].append(*(thisFile->getMapLabelTable(i)), true);//produce an error on conflict
                }
            }
            if (!discardUnusedLabels)//Matt insisted on making this an error, too
            {
                if (myXML.getNumberOfDimensions() != 2) throw AlgorithmException("only 2D cifti are supported");
                if (myXML.getMappingType(1 - myDir) != CiftiMappingType::LABELS) throw AlgorithmException("cannot replace structure with label data in non-label cifti file");
                const CiftiLabelsMap& myLabelsMap = myXML.getLabelsMap(1 - myDir);
                if (myLabelsMap.getLength() != int64_t(testTables.size())) throw AlgorithmException("input files have the wrong number of maps for this cifti file");
                try
                {
                    for (int i = 0; i < int(testTables.size()); ++i)
                    {
                        testTables[i].append(*(myLabelsMap.getMapLabelTable(i)), true);
                    }
                } catch (CaretException& e) {
                    throw AlgorithmException("labels in the cifti file conflict with labels in other input files, but all structures are being replaced, so using either -label-collision or -discard-unused-labels will resolve this error");
                }
            }
        }
    }
    vector<DataSourceInfo> grayordOrder;
    {
        //most predictable is probably to do -cifti-create-label order and ignore whatever the input cifti order is
        auto iter = surfSources.find(StructureEnum::CORTEX_LEFT);
        if (iter != surfSources.end())
        {
            grayordOrder.push_back(iter->second);
            surfSources.erase(iter);
        }
        iter = surfSources.find(StructureEnum::CORTEX_RIGHT);
        if (iter != surfSources.end())
        {
            grayordOrder.push_back(iter->second);
            surfSources.erase(iter);
        }
        //cerebellum next, because this is what -cifti-create-label does
        iter = surfSources.find(StructureEnum::CEREBELLUM);
        if (iter != surfSources.end())
        {
            grayordOrder.push_back(iter->second);
            surfSources.erase(iter);
        }
        //volume after main surfaces
        if (volumeAllOpt->m_present)
        {
            grayordOrder.push_back(DataSourceInfo(VOLUME_ALL, -1));
        }
        for (auto vIter : volSources) //this is already checked as mutually exclusive with -volume-all
        {
            grayordOrder.push_back(vIter.second);
        }
        volSources.clear();
        //any other surfaces last?  -cifti-create-label doesn't accept these
        for (auto sIter : surfSources)
        {
            grayordOrder.push_back(sIter.second);
        }
        surfSources.clear();
        CaretAssert(volSources.empty() && surfSources.empty());
    }
    vector<DataSourceInfo> useOrder;
    switch (conflictLogic)
    {
        case ERROR://do scalar-type replace with legacy order, in case it makes a difference
        case LEGACY:
            useOrder = legacyOrder;
            break;
        case LEFT_SURFACE_FIRST:
            useOrder = vector<DataSourceInfo>(grayordOrder.rbegin(), grayordOrder.rend());//reverse the order so high priority goes last
            break;
    }
    for (auto iter : useOrder)
    {
        switch (iter.type)
        {
            case LABEL:
            {
                AString structName = labelInst[iter.index]->getString(1);
                bool ok = false;
                StructureEnum::Enum myStruct = StructureEnum::fromName(structName, &ok);
                if (!ok)
                {
                    throw AlgorithmException("unrecognized structure name");
                }
                LabelFile* labelIn = labelInst[iter.index]->getLabel(2);
                AlgorithmCiftiReplaceStructure(NULL, &myCifti, myDir, myStruct, labelIn, discardUnusedLabels, errorOnLabelConflict);
                break;
            }
            case METRIC:
            {
                AString structName = metricInst[iter.index]->getString(1);
                bool ok = false;
                StructureEnum::Enum myStruct = StructureEnum::fromName(structName, &ok);
                if (!ok)
                {
                    throw AlgorithmException("unrecognized structure name");
                }
                MetricFile* metricIn = metricInst[iter.index]->getMetric(2);
                AlgorithmCiftiReplaceStructure(NULL, &myCifti, myDir, myStruct, metricIn);
                break;
            }
            case VOLUME:
            {//can't happen with -volume-all
                AString structName = volumeInst[iter.index]->getString(1);
                bool ok = false;
                StructureEnum::Enum myStruct = StructureEnum::fromName(structName, &ok);
                if (!ok)
                {
                    throw AlgorithmException("unrecognized structure name");
                }
                VolumeFile* volIn = volumeInst[iter.index]->getVolume(2);
                bool fromCropVol = volumeInst[iter.index]->getOptionalParameter(3)->m_present;
                AlgorithmCiftiReplaceStructure(NULL, &myCifti, myDir, myStruct, volIn, fromCropVol, discardUnusedLabels, errorOnLabelConflict);
                break;
            }
            case VOLUME_ALL:
            {//can only exist once
                VolumeFile* volIn = volumeAllOpt->getVolume(1);
                bool fromCropVol = volumeAllOpt->getOptionalParameter(2)->m_present;
                AlgorithmCiftiReplaceStructure(NULL, &myCifti, myDir, volIn, fromCropVol, discardUnusedLabels, errorOnLabelConflict);
                break;
            }
        }
    }
    myCifti.writeFile(ciftiName);//and write the modified file
}

AlgorithmCiftiReplaceStructure::AlgorithmCiftiReplaceStructure(ProgressObject* myProgObj, CiftiFile* ciftiInOut, const int myDir,
                                                               const StructureEnum::Enum myStruct, const MetricFile* metricIn) : AbstractAlgorithm(myProgObj)
{
    LevelProgress myProgress(myProgObj);
    const CiftiXML& myXML = ciftiInOut->getCiftiXML();
    if (myXML.getNumberOfDimensions() != 2) throw AlgorithmException("replace structure only supported on 2D cifti");
    if (myXML.getMappingType(myDir) != CiftiMappingType::BRAIN_MODELS) throw AlgorithmException("specified dimension must contain brain models");
    if (myXML.getMappingType(1 - myDir) == CiftiMappingType::LABELS) CaretLogWarning("replacing data in label-type cifti file with a metric file!");//previously not an error, so just warn, I guess
    const CiftiBrainModelsMap& myDenseMap = myXML.getBrainModelsMap(myDir);
    vector<CiftiBrainModelsMap::SurfaceMap> myMap;
    int rowSize = ciftiInOut->getNumberOfColumns(), colSize = ciftiInOut->getNumberOfRows();
    if (myDir == CiftiXML::ALONG_COLUMN)
    {
        myMap = myDenseMap.getSurfaceMap(myStruct);
        int64_t numNodes = myDenseMap.getSurfaceNumberOfNodes(myStruct);
        if (metricIn->getNumberOfNodes() != numNodes)
        {
            throw AlgorithmException("input metric has the wrong number of vertices");
        }
        if (metricIn->getNumberOfColumns() != rowSize)
        {
            throw AlgorithmException("input metric has the wrong number of columns");
        }
        int mapSize = (int)myMap.size();
        CaretArray<float> rowScratch(rowSize);
        for (int i = 0; i < mapSize; ++i)
        {
            for (int j = 0; j < rowSize; ++j)
            {
                rowScratch[j] = metricIn->getValue(myMap[i].m_surfaceNode, j);
            }
            ciftiInOut->setRow(rowScratch, myMap[i].m_ciftiIndex);
        }
    } else {
        if (myDir != CiftiXML::ALONG_ROW) throw AlgorithmException("unsupported cifti direction");
        myMap = myDenseMap.getSurfaceMap(myStruct);
        int64_t numNodes = myDenseMap.getSurfaceNumberOfNodes(myStruct);
        if (metricIn->getNumberOfNodes() != numNodes)
        {
            throw AlgorithmException("input metric has the wrong number of vertices");
        }
        if (metricIn->getNumberOfColumns() != colSize)
        {
            throw AlgorithmException("input metric has the wrong number of columns");
        }
        int mapSize = (int)myMap.size();
        CaretArray<float> rowScratch(rowSize);
        for (int i = 0; i < colSize; ++i)
        {
            ciftiInOut->getRow(rowScratch, i, true);
            for (int j = 0; j < mapSize; ++j)
            {
                rowScratch[myMap[j].m_ciftiIndex] = metricIn->getValue(myMap[j].m_surfaceNode, i);
            }
            ciftiInOut->setRow(rowScratch, i);
        }
    }
}

AlgorithmCiftiReplaceStructure::AlgorithmCiftiReplaceStructure(ProgressObject* myProgObj, CiftiFile* ciftiInOut, const int myDir,
                                                               const StructureEnum::Enum myStruct, const LabelFile* labelIn,
                                                               const bool discardUnusedLabels, const bool errorOnLabelConflict) : AbstractAlgorithm(myProgObj)
{
    LevelProgress myProgress(myProgObj);
    const CiftiXML& myXML = ciftiInOut->getCiftiXML();
    if (myXML.getNumberOfDimensions() != 2) throw AlgorithmException("replace structure only supported on 2D cifti");
    if (myXML.getMappingType(myDir) != CiftiMappingType::BRAIN_MODELS) throw AlgorithmException("specified dimension must contain brain models");
    const CiftiBrainModelsMap& myDenseMap = myXML.getBrainModelsMap(myDir);
    ciftiInOut->convertToInMemory();//so that writing it can change the header
    vector<CiftiBrainModelsMap::SurfaceMap> myMap;
    int64_t rowSize = ciftiInOut->getNumberOfColumns(), colSize = ciftiInOut->getNumberOfRows();
    if (myDir == CiftiXML::ALONG_COLUMN)
    {
        if (myXML.getMappingType(CiftiXML::ALONG_ROW) != CiftiMappingType::LABELS) throw AlgorithmException("label replace structure requested on non-label cifti");
        const CiftiLabelsMap& myLabelsMap = myXML.getLabelsMap(CiftiXML::ALONG_ROW);
        myMap = myDenseMap.getSurfaceMap(myStruct);
        int64_t numNodes = myDenseMap.getSurfaceNumberOfNodes(myStruct);
        if (labelIn->getNumberOfNodes() != numNodes)
        {
            throw AlgorithmException("input label file has the wrong number of vertices");
        }
        if (labelIn->getNumberOfColumns() != rowSize)
        {
            throw AlgorithmException("input label file has the wrong number of columns");
        }
        int64_t mapSize = (int64_t)myMap.size();
        CaretArray<float> rowScratch(rowSize, 0.0f);//initialize to the default "unused" value in case we get short reads due to unallocated on-disk cifti
        vector<set<int32_t> > usedArray(rowSize);
        vector<map<int32_t, int32_t> > remapArray(rowSize);
        for (int64_t j = 0; j < rowSize; ++j)
        {
            GiftiLabelTable myTable = *(labelIn->getLabelTable());//we remap the old label table values so that the new label table keys are unmolested
            remapArray[j] = myTable.append(*(myLabelsMap.getMapLabelTable(j)), errorOnLabelConflict);
            *(myLabelsMap.getMapLabelTable(j)) = myTable;
        }
        set<int64_t> writeRows;
        for (int64_t i = 0; i < mapSize; ++i)
        {
            writeRows.insert(myMap[i].m_ciftiIndex);
        }
        for (int64_t i = 0; i < colSize; ++i)//we have to remap all old data since we are changing he keys of the old label table
        {
            if (writeRows.find(i) == writeRows.end())
            {
                ciftiInOut->getRow(rowScratch, i, true);
                bool rowChanged = false;
                for (int64_t j = 0; j < rowSize; ++j)
                {
                    int32_t tempKey = (int32_t)floor(rowScratch[j] + 0.5f);
                    map<int32_t, int32_t>::iterator iter = remapArray[j].find(tempKey);
                    if (iter != remapArray[j].end())
                    {
                        if (iter->first != iter->second)
                        {
                            rowChanged = true;
                            rowScratch[j] = iter->second;//remap the key if its value changes
                            usedArray[j].insert(iter->second);
                        } else {
                            usedArray[j].insert(tempKey);
                        }
                    } else {
                        rowChanged = true;
                        int32_t unusedLabel = myLabelsMap.getMapLabelTable(j)->getUnassignedLabelKey();
                        rowScratch[j] = unusedLabel;
                        usedArray[j].insert(unusedLabel);
                    }
                }
                if (rowChanged) ciftiInOut->setRow(rowScratch, i);
            }
        }
        for (int64_t i = 0; i < mapSize; ++i)//set the new rows
        {
            for (int64_t j = 0; j < rowSize; ++j)
            {
                int32_t tempKey = labelIn->getLabelKey(myMap[i].m_surfaceNode, j);
                usedArray[j].insert(tempKey);
                rowScratch[j] = tempKey;
            }
            ciftiInOut->setRow(rowScratch, myMap[i].m_ciftiIndex);
        }
        if (discardUnusedLabels)
        {
            for (int64_t i = 0; i < rowSize; ++i)//delete unused labels
            {
                myLabelsMap.getMapLabelTable(i)->deleteUnusedLabels(usedArray[i]);
            }
        }
    } else {
        if (myDir != CiftiXML::ALONG_ROW) throw AlgorithmException("unsupported cifti direction");
        if (myXML.getMappingType(CiftiXML::ALONG_COLUMN) != CiftiMappingType::LABELS) throw AlgorithmException("label replace structure requested on non-label cifti");
        const CiftiLabelsMap& myLabelsMap = myXML.getLabelsMap(CiftiXML::ALONG_COLUMN);
        myMap = myDenseMap.getSurfaceMap(myStruct);
        {
            throw AlgorithmException("structure not found in specified dimension");
        }
        int64_t numNodes = myDenseMap.getSurfaceNumberOfNodes(myStruct);
        if (labelIn->getNumberOfNodes() != numNodes)
        {
            throw AlgorithmException("input label file has the wrong number of vertices");
        }
        if (labelIn->getNumberOfColumns() != colSize)
        {
            throw AlgorithmException("input label file has the wrong number of columns");
        }
        int64_t mapSize = (int64_t)myMap.size();
        CaretArray<float> rowScratch(rowSize, 0.0f);
        set<int64_t> writeCols;
        for (int64_t i = 0; i < mapSize; ++i)
        {
            writeCols.insert(myMap[i].m_ciftiIndex);
        }
        for (int64_t i = 0; i < colSize; ++i)
        {
            GiftiLabelTable myTable = *(labelIn->getLabelTable());//we remap the old label table values so that the new label table keys are unmolested
            map<int32_t, int32_t> remap = myTable.append(*(myLabelsMap.getMapLabelTable(i)), errorOnLabelConflict);
            *(myLabelsMap.getMapLabelTable(i)) = myTable;
            ciftiInOut->getRow(rowScratch, i, true);
            set<int32_t> used;
            int32_t unusedKey = myTable.getUnassignedLabelKey();
            for(int64_t j = 0; j < rowSize; ++j)
            {
                if (writeCols.find(j) == writeCols.end())
                {
                    int32_t tempKey = (int32_t)floor(rowScratch[j] + 0.5f);
                    map<int32_t, int32_t>::iterator iter = remap.find(tempKey);
                    if (iter != remap.end())
                    {
                        rowScratch[j] = iter->second;
                        used.insert(iter->second);
                    } else {
                        rowScratch[j] = unusedKey;
                        used.insert(unusedKey);
                    }
                }
            }
            for (int64_t j = 0; j < mapSize; ++j)
            {
                int32_t tempKey = labelIn->getLabelKey(myMap[j].m_surfaceNode, i);
                rowScratch[myMap[j].m_ciftiIndex] = tempKey;
                used.insert(tempKey);
            }
            if (discardUnusedLabels)
            {
                myLabelsMap.getMapLabelTable(i)->deleteUnusedLabels(used);
            }
            ciftiInOut->setRow(rowScratch, i);
        }
    }
}

AlgorithmCiftiReplaceStructure::AlgorithmCiftiReplaceStructure(ProgressObject* myProgObj, CiftiFile* ciftiInOut, const int myDir,
                                                               const StructureEnum::Enum myStruct, const VolumeFile* volIn, const bool fromCropped,
                                                               const bool discardUnusedLabels, const bool errorOnLabelConflict) : AbstractAlgorithm(myProgObj)
{
    const CiftiXML& myXML = ciftiInOut->getCiftiXML();
    if (myDir != CiftiXML::ALONG_ROW && myDir != CiftiXML::ALONG_COLUMN) throw AlgorithmException("direction not supported in cifti replace structure");
    if (myXML.getNumberOfDimensions() != 2) throw AlgorithmException("replace structure only supported on 2D cifti");
    LevelProgress myProgress(myProgObj);
    int64_t rowSize = ciftiInOut->getNumberOfColumns(), colSize = ciftiInOut->getNumberOfRows();
    if (myXML.getMappingType(myDir) != CiftiMappingType::BRAIN_MODELS) throw AlgorithmException("specified direction does not contain brain models");
    const CiftiBrainModelsMap& myBrainMap = myXML.getBrainModelsMap(myDir);
    const int64_t* myDims = myBrainMap.getVolumeSpace().getDims();
    vector<vector<float> > mySform = myBrainMap.getVolumeSpace().getSform();
    vector<CiftiBrainModelsMap::VolumeMap> myMap = myBrainMap.getVolumeStructureMap(myStruct);
    int64_t numVoxels = (int64_t)myMap.size();
    int64_t offset[3];
    vector<int64_t> newdims;
    if (fromCropped)
    {
        newdims.resize(3);
        AlgorithmCiftiSeparate::getCroppedVolSpace(ciftiInOut, myDir, myStruct, newdims.data(), mySform, offset);
    } else {
        newdims.push_back(myDims[0]);
        newdims.push_back(myDims[1]);
        newdims.push_back(myDims[2]);
        offset[0] = 0;
        offset[1] = 0;
        offset[2] = 0;
    }
    if (!volIn->matchesVolumeSpace(newdims.data(), mySform))
    {
        throw AlgorithmException("input volume doesn't match volume space and dimensions in CIFTI");
    }
    CaretArray<float> rowScratch(rowSize, 0.0f);//initialize with default unused label in case dlabel and short reads due to uninitialized on-disk cifti
    vector<int64_t> volDims;
    volIn->getDimensions(volDims);
    if (myDir == CiftiXML::ALONG_COLUMN)
    {
        if (volDims[3] != rowSize)
        {
            throw AlgorithmException("volume has the wrong number of subvolumes");
        }
        if (myXML.getMappingType(CiftiXML::ALONG_ROW) == CiftiMappingType::LABELS)
        {
            if (volIn->getType() != SubvolumeAttributes::LABEL) throw AlgorithmException("replace structure called on cifti label file with non-label volume");
            ciftiInOut->convertToInMemory();//so that writing can change the XML
            const CiftiLabelsMap& myLabelsMap = myXML.getLabelsMap(CiftiXML::ALONG_ROW);
            vector<map<int32_t, int32_t> > remapArray(rowSize);
            vector<set<int32_t> > usedArray(rowSize);
            for (int64_t i = 0; i < rowSize; ++i)
            {
                GiftiLabelTable myTable = *(volIn->getMapLabelTable(i));//we remap the old label table values so that the new label table keys are unmolested
                remapArray[i] = myTable.append(*(myLabelsMap.getMapLabelTable(i)), errorOnLabelConflict);
                *(myLabelsMap.getMapLabelTable(i)) = myTable;
            }
            set<int64_t> writeRows;
            for (int64_t i = 0; i < numVoxels; ++i)
            {
                writeRows.insert(myMap[i].m_ciftiIndex);
            }
            for (int64_t i = 0; i < colSize; ++i)
            {
                if (writeRows.find(i) == writeRows.end())
                {
                    ciftiInOut->getRow(rowScratch, i, true);
                    bool rowChanged = false;
                    for (int64_t j = 0; j < rowSize; ++j)
                    {
                        int32_t tempKey = (int32_t)floor(rowScratch[j] + 0.5f);
                        map<int32_t, int32_t>::iterator iter = remapArray[j].find(tempKey);
                        if (iter != remapArray[j].end())
                        {
                            if (iter->first != iter->second)
                            {
                                rowChanged = true;
                                rowScratch[j] = iter->second;//remap the key if its value changes
                                usedArray[j].insert(iter->second);
                            } else {
                                usedArray[j].insert(tempKey);
                            }
                        } else {
                            rowChanged = true;
                            int32_t unusedLabel = myLabelsMap.getMapLabelTable(j)->getUnassignedLabelKey();
                            rowScratch[j] = unusedLabel;
                            usedArray[j].insert(unusedLabel);
                        }
                    }
                    if (rowChanged) ciftiInOut->setRow(rowScratch, i);
                }
            }
            for (int64_t i = 0; i < numVoxels; ++i)
            {
                int64_t thisvoxel[3] = { myMap[i].m_ijk[0] - offset[0], myMap[i].m_ijk[1] - offset[1], myMap[i].m_ijk[2] - offset[2] };
                for (int j = 0; j < rowSize; ++j)
                {
                    int32_t tempKey = (int32_t)floor(volIn->getValue(thisvoxel, j) + 0.5f);
                    usedArray[j].insert(tempKey);
                    rowScratch[j] = tempKey;
                }
                ciftiInOut->setRow(rowScratch, myMap[i].m_ciftiIndex);
            }
            if (discardUnusedLabels)
            {
                for (int64_t i = 0; i < rowSize; ++i)//delete unused labels
                {
                    myLabelsMap.getMapLabelTable(i)->deleteUnusedLabels(usedArray[i]);
                }
            }
        } else {
            for (int64_t i = 0; i < numVoxels; ++i)
            {
                int64_t thisvoxel[3] = { myMap[i].m_ijk[0] - offset[0], myMap[i].m_ijk[1] - offset[1], myMap[i].m_ijk[2] - offset[2] };
                for (int j = 0; j < rowSize; ++j)
                {
                    rowScratch[j] = volIn->getValue(thisvoxel, j);
                }
                ciftiInOut->setRow(rowScratch, myMap[i].m_ciftiIndex);
            }
        }
    } else {
        if (volDims[3] != colSize)
        {
            throw AlgorithmException("volume has the wrong number of subvolumes");
        }
        if (myXML.getMappingType(CiftiXML::ALONG_COLUMN) == CiftiMappingType::LABELS)
        {
            if (volIn->getType() != SubvolumeAttributes::LABEL) throw AlgorithmException("replace structure called on cifti label file with non-label volume");
            ciftiInOut->convertToInMemory();
            const CiftiLabelsMap& myLabelsMap = myXML.getLabelsMap(CiftiXML::ALONG_COLUMN);
            set<int64_t> writeCols;
            for (int64_t i = 0; i < numVoxels; ++i)
            {
                writeCols.insert(myMap[i].m_ciftiIndex);
            }
            for (int64_t i = 0; i < colSize; ++i)
            {
                GiftiLabelTable myTable = *(volIn->getMapLabelTable(i));//we remap the old label table values so that the new label table keys are unmolested
                map<int32_t, int32_t> remap = myTable.append(*(myLabelsMap.getMapLabelTable(i)), errorOnLabelConflict);
                *(myLabelsMap.getMapLabelTable(i)) = myTable;
                ciftiInOut->getRow(rowScratch, i, true);
                set<int32_t> used;
                int32_t unusedKey = myTable.getUnassignedLabelKey();
                for(int64_t j = 0; j < rowSize; ++j)
                {
                    if (writeCols.find(j) == writeCols.end())
                    {
                        int32_t tempKey = (int32_t)floor(rowScratch[j] + 0.5f);
                        map<int32_t, int32_t>::iterator iter = remap.find(tempKey);
                        if (iter != remap.end())
                        {
                            rowScratch[j] = iter->second;
                            used.insert(iter->second);
                        } else {
                            rowScratch[j] = unusedKey;
                            used.insert(unusedKey);
                        }
                    }
                }
                for (int64_t j = 0; j < numVoxels; ++j)
                {
                    int64_t thisvoxel[3] = { myMap[j].m_ijk[0] - offset[0], myMap[j].m_ijk[1] - offset[1], myMap[j].m_ijk[2] - offset[2] };
                    int32_t tempKey = (int32_t)floor(volIn->getValue(thisvoxel, i) + 0.5f);
                    rowScratch[myMap[j].m_ciftiIndex] = tempKey;
                    used.insert(tempKey);
                }
                if (discardUnusedLabels)
                {
                    myLabelsMap.getMapLabelTable(i)->deleteUnusedLabels(used);
                }
                ciftiInOut->setRow(rowScratch, i);
            }
        } else {
            for (int64_t i = 0; i < colSize; ++i)
            {
                ciftiInOut->getRow(rowScratch, i, true);//the on-disk cifti file may not have been allocated yet, so short reads are okay
                for (int64_t j = 0; j < numVoxels; ++j)
                {
                    int64_t thisvoxel[3] = { myMap[j].m_ijk[0] - offset[0], myMap[j].m_ijk[1] - offset[1], myMap[j].m_ijk[2] - offset[2] };
                    rowScratch[myMap[j].m_ciftiIndex] = volIn->getValue(thisvoxel, i);
                }
                ciftiInOut->setRow(rowScratch, i);
            }
        }
    }
}

AlgorithmCiftiReplaceStructure::AlgorithmCiftiReplaceStructure(ProgressObject* myProgObj, CiftiFile* ciftiInOut, const int myDir,
                                                               const VolumeFile* volIn, const bool fromCropped,
                                                               const bool discardUnusedLabels, const bool errorOnLabelConflict): AbstractAlgorithm(myProgObj)
{
    const CiftiXML& myXML = ciftiInOut->getCiftiXML();
    if (myXML.getNumberOfDimensions() != 2) throw AlgorithmException("replace structure only supported on 2D cifti");
    LevelProgress myProgress(myProgObj);
    if (myDir != CiftiXML::ALONG_ROW && myDir != CiftiXML::ALONG_COLUMN) throw AlgorithmException("direction not supported in cifti replace structure");
    if (myXML.getMappingType(myDir) != CiftiMappingType::BRAIN_MODELS) throw AlgorithmException("specified direction does not contain brain models");
    const CiftiBrainModelsMap& myBrainMap = myXML.getBrainModelsMap(myDir);
    const int64_t* myDims = myBrainMap.getVolumeSpace().getDims();
    vector<vector<float> > mySform = myBrainMap.getVolumeSpace().getSform();
    vector<CiftiBrainModelsMap::VolumeMap> myMap = myBrainMap.getFullVolumeMap();
    int64_t numVoxels = (int64_t)myMap.size();
    int64_t rowSize = ciftiInOut->getNumberOfColumns(), colSize = ciftiInOut->getNumberOfRows();
    vector<int64_t> newdims;
    int64_t offset[3];
    if (fromCropped)
    {
        newdims.resize(3);
        AlgorithmCiftiSeparate::getCroppedVolSpaceAll(ciftiInOut, myDir, newdims.data(), mySform, offset);
    } else {
        newdims.push_back(myDims[0]);
        newdims.push_back(myDims[1]);
        newdims.push_back(myDims[2]);
        offset[0] = 0;
        offset[1] = 0;
        offset[2] = 0;
    }
    if (!volIn->matchesVolumeSpace(newdims.data(), mySform))
    {
        throw AlgorithmException("input volume doesn't match volume space and dimensions in CIFTI");
    }
    CaretArray<float> rowScratch(rowSize, 0.0f);//initialize with default unused label in case dlabel and short reads due to uninitialized on-disk cifti
    vector<int64_t> volDims;
    volIn->getDimensions(volDims);
    if (myDir == CiftiXML::ALONG_COLUMN)
    {
        if (volDims[3] != rowSize)
        {
            throw AlgorithmException("volume has the wrong number of subvolumes");
        }
        if (myXML.getMappingType(CiftiXML::ALONG_ROW) == CiftiMappingType::LABELS)
        {
            const CiftiLabelsMap& myLabelsMap = myXML.getLabelsMap(CiftiXML::ALONG_ROW);
            if (volIn->getType() != SubvolumeAttributes::LABEL) throw AlgorithmException("replace structure called on cifti label file with non-label volume");
            ciftiInOut->convertToInMemory();
            vector<map<int32_t, int32_t> > remapArray(rowSize);
            vector<set<int32_t> > usedArray(rowSize);
            for (int64_t i = 0; i < rowSize; ++i)
            {
                GiftiLabelTable myTable = *(volIn->getMapLabelTable(i));//we remap the old label table values so that the new label table keys are unmolested
                remapArray[i] = myTable.append(*(myLabelsMap.getMapLabelTable(i)), errorOnLabelConflict);
                *(myLabelsMap.getMapLabelTable(i)) = myTable;
            }
            set<int64_t> writeRows;
            for (int64_t i = 0; i < numVoxels; ++i)
            {
                writeRows.insert(myMap[i].m_ciftiIndex);
            }
            for (int64_t i = 0; i < colSize; ++i)
            {
                if (writeRows.find(i) == writeRows.end())
                {
                    ciftiInOut->getRow(rowScratch, i, true);
                    bool rowChanged = false;
                    for (int64_t j = 0; j < rowSize; ++j)
                    {
                        int32_t tempKey = (int32_t)floor(rowScratch[j] + 0.5f);
                        map<int32_t, int32_t>::iterator iter = remapArray[j].find(tempKey);
                        if (iter != remapArray[j].end())
                        {
                            if (iter->first != iter->second)
                            {
                                rowChanged = true;
                                rowScratch[j] = iter->second;//remap the key if its value changes
                                usedArray[j].insert(iter->second);
                            } else {
                                usedArray[j].insert(tempKey);
                            }
                        } else {
                            rowChanged = true;
                            int32_t unusedLabel = myLabelsMap.getMapLabelTable(j)->getUnassignedLabelKey();
                            rowScratch[j] = unusedLabel;
                            usedArray[j].insert(unusedLabel);
                        }
                    }
                    if (rowChanged) ciftiInOut->setRow(rowScratch, i);
                }
            }
            for (int64_t i = 0; i < numVoxels; ++i)
            {
                int64_t thisvoxel[3] = { myMap[i].m_ijk[0] - offset[0], myMap[i].m_ijk[1] - offset[1], myMap[i].m_ijk[2] - offset[2] };
                for (int j = 0; j < rowSize; ++j)
                {
                    int32_t tempKey = (int32_t)floor(volIn->getValue(thisvoxel, j) + 0.5f);
                    usedArray[j].insert(tempKey);
                    rowScratch[j] = tempKey;
                }
                ciftiInOut->setRow(rowScratch, myMap[i].m_ciftiIndex);
            }
            if (discardUnusedLabels)
            {
                for (int64_t i = 0; i < rowSize; ++i)//delete unused labels
                {
                    myLabelsMap.getMapLabelTable(i)->deleteUnusedLabels(usedArray[i]);
                }
            }
        } else {
            for (int64_t i = 0; i < numVoxels; ++i)
            {
                int64_t thisvoxel[3] = { myMap[i].m_ijk[0] - offset[0], myMap[i].m_ijk[1] - offset[1], myMap[i].m_ijk[2] - offset[2] };
                for (int j = 0; j < rowSize; ++j)
                {
                    rowScratch[j] = volIn->getValue(thisvoxel, j);
                }
                ciftiInOut->setRow(rowScratch, myMap[i].m_ciftiIndex);
            }
        }
    } else {
        if (volDims[3] != colSize)
        {
            throw AlgorithmException("volume has the wrong number of subvolumes");
        }
        if (myXML.getMappingType(CiftiXML::ALONG_COLUMN) == CiftiMappingType::LABELS)
        {
            const CiftiLabelsMap& myLabelsMap = myXML.getLabelsMap(CiftiXML::ALONG_COLUMN);
            if (volIn->getType() != SubvolumeAttributes::LABEL) throw AlgorithmException("replace structure called on cifti label file with non-label volume");
            ciftiInOut->convertToInMemory();
            set<int64_t> writeCols;
            for (int64_t i = 0; i < numVoxels; ++i)
            {
                writeCols.insert(myMap[i].m_ciftiIndex);
            }
            for (int64_t i = 0; i < colSize; ++i)
            {
                GiftiLabelTable myTable = *(volIn->getMapLabelTable(i));//we remap the old label table values so that the new label table keys are unmolested
                map<int32_t, int32_t> remap = myTable.append(*(myLabelsMap.getMapLabelTable(i)), errorOnLabelConflict);
                *(myLabelsMap.getMapLabelTable(i)) = myTable;
                ciftiInOut->getRow(rowScratch, i, true);
                set<int32_t> used;
                int32_t unusedKey = myTable.getUnassignedLabelKey();
                for(int64_t j = 0; j < rowSize; ++j)
                {
                    if (writeCols.find(j) == writeCols.end())
                    {
                        int32_t tempKey = (int32_t)floor(rowScratch[j] + 0.5f);
                        map<int32_t, int32_t>::iterator iter = remap.find(tempKey);
                        if (iter != remap.end())
                        {
                            rowScratch[j] = iter->second;
                            used.insert(iter->second);
                        } else {
                            rowScratch[j] = unusedKey;
                            used.insert(unusedKey);
                        }
                    }
                }
                for (int64_t j = 0; j < numVoxels; ++j)
                {
                    int64_t thisvoxel[3] = { myMap[j].m_ijk[0] - offset[0], myMap[j].m_ijk[1] - offset[1], myMap[j].m_ijk[2] - offset[2] };
                    int32_t tempKey = (int32_t)floor(volIn->getValue(thisvoxel, i) + 0.5f);
                    rowScratch[myMap[j].m_ciftiIndex] = tempKey;
                    used.insert(tempKey);
                }
                if (discardUnusedLabels)
                {
                    myLabelsMap.getMapLabelTable(i)->deleteUnusedLabels(used);
                }
                ciftiInOut->setRow(rowScratch, i);
            }
        } else {
            for (int64_t i = 0; i < colSize; ++i)
            {
                ciftiInOut->getRow(rowScratch, i, true);//the on-disk cifti file may not have been allocated yet, so short reads are okay
                for (int64_t j = 0; j < numVoxels; ++j)
                {
                    int64_t thisvoxel[3] = { myMap[j].m_ijk[0] - offset[0], myMap[j].m_ijk[1] - offset[1], myMap[j].m_ijk[2] - offset[2] };
                    rowScratch[myMap[j].m_ciftiIndex] = volIn->getValue(thisvoxel, i);
                }
                ciftiInOut->setRow(rowScratch, i);
            }
        }
    }
}

float AlgorithmCiftiReplaceStructure::getAlgorithmInternalWeight()
{
    return 1.0f;//override this if needed, if the progress bar isn't smooth
}

float AlgorithmCiftiReplaceStructure::getSubAlgorithmWeight()
{
    return 0.0f;
}
