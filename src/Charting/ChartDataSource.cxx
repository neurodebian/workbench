
/*LICENSE_START*/
/*
 * Copyright 2014 Washington University,
 * All rights reserved.
 *
 * Connectome DB and Connectome Workbench are part of the integrated Connectome
 * Informatics Platform.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *    * Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *    * Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *    * Neither the names of Washington University nor the
 *      names of its contributors may be used to endorse or promote products
 *      derived from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR  SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */
/*LICENSE_END*/

#include <QXmlStreamReader>
#include <QXmlStreamWriter>

#define __CHART_DATA_SOURCE_DECLARE__
#include "ChartDataSource.h"
#undef __CHART_DATA_SOURCE_DECLARE__

#include "CaretAssert.h"
#include "CaretLogger.h"
#include "SceneClass.h"
#include "SceneClassArray.h"
#include "SceneClassAssistant.h"

using namespace caret;



/**
 * \class caret::ChartDataSource
 * \brief Contains source of data that is displayed in a chart.
 * \ingroup Charting
 */

/**
 * Constructor.
 */
ChartDataSource::ChartDataSource()
: CaretObject()
{
    initializeMembersChartDataSource();
}

/**
 * Destructor.
 */
ChartDataSource::~ChartDataSource()
{
    delete m_sceneAssistant;
}

/**
 * Initialize members of a new instance.
 */
void
ChartDataSource::initializeMembersChartDataSource()
{
    m_dataSourceMode = ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_INVALID;
    m_nodeIndex = -1;
    m_voxelIJK[0] = -1;
    m_voxelIJK[1] = -1;
    m_voxelIJK[2] = -1;
    
    m_sceneAssistant = new SceneClassAssistant();
    m_sceneAssistant->add<ChartDataSourceModeEnum, ChartDataSourceModeEnum::Enum>("m_dataSourceMode",
                                                                                  &m_dataSourceMode);
    m_sceneAssistant->add("m_nodeIndex",
                          &m_nodeIndex);
    m_sceneAssistant->add("m_surfaceNumberOfNodes",
                          &m_surfaceNumberOfNodes);
    m_sceneAssistant->add("m_surfaceStructureName",
                          &m_surfaceStructureName);
    m_sceneAssistant->addArray("m_voxelIJK",
                               m_voxelIJK, 3, -1);
}

/**
 * Copy constructor.
 * @param obj
 *    Object that is copied.
 */
ChartDataSource::ChartDataSource(const ChartDataSource& obj)
: CaretObject(obj)
{
    initializeMembersChartDataSource();
    
    this->copyHelperChartDataSource(obj);
}

/**
 * Assignment operator.
 * @param obj
 *    Data copied from obj to this.
 * @return
 *    Reference to this object.
 */
ChartDataSource&
ChartDataSource::operator=(const ChartDataSource& obj)
{
    if (this != &obj) {
        CaretObject::operator=(obj);
        this->copyHelperChartDataSource(obj);
    }
    return *this;
}

/**
 * Helps with copying an object of this type.
 * @param obj
 *    Object that is copied.
 */
void
ChartDataSource::copyHelperChartDataSource(const ChartDataSource& obj)
{
    m_dataSourceMode       = obj.m_dataSourceMode;
    m_chartableFileName    = obj.m_chartableFileName;
    m_nodeIndex            = obj.m_nodeIndex;
    m_surfaceNumberOfNodes = obj.m_surfaceNumberOfNodes;
    m_surfaceStructureName = obj.m_surfaceStructureName;
    m_nodeIndicesAverage   = obj.m_nodeIndicesAverage;
    m_voxelIJK[0]          = obj.m_voxelIJK[0];
    m_voxelIJK[1]          = obj.m_voxelIJK[1];
    m_voxelIJK[2]          = obj.m_voxelIJK[2];
}

/**
 * @return Name of the chartable file.
 */
AString
ChartDataSource::getChartableFileName() const
{
    return m_chartableFileName;
}


/**
 * Setup for a surface node source.
 *
 * @param chartableFileName
 *    Name of the chartable file.
 * @param surfaceStructureName
 *    Name of surface structure.
 * @param surfaceNumberOfNodes
 *    Number of nodes in the surface.
 * @param nodeIndex
 *    Index of the surface node.
 */
void
ChartDataSource::setSurfaceNode(const AString& chartableFileName,
                                const AString& surfaceStructureName,
                                const int32_t surfaceNumberOfNodes,
                                const int32_t nodeIndex)
{
    CaretAssert(nodeIndex >= 0);
    
    m_dataSourceMode = ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_SURFACE_NODE_INDEX;
    m_chartableFileName = chartableFileName;
    m_surfaceNumberOfNodes = surfaceNumberOfNodes;
    m_surfaceStructureName = surfaceStructureName;
    m_nodeIndex = nodeIndex;
}

/**
 * @return Mode indicating source of the data.
 */
ChartDataSourceModeEnum::Enum
ChartDataSource::getDataSourceMode() const
{
    return m_dataSourceMode;
}

/**
 * Get the surface node data source.
 *
 * @param surfaceStructureName
 *    Name of surface structure.
 * @param surfaceNumberOfNodes
 *    Number of nodes in the surface.
 * @param nodeIndex
 *    Index of the surface node.
 */
void
ChartDataSource::getSurfaceNode(AString& surfaceStructureName,
                                int32_t& surfaceNumberOfNodes,
                                int32_t& nodeIndex) const
{
    surfaceStructureName = m_surfaceStructureName;
    surfaceNumberOfNodes = m_surfaceNumberOfNodes;
    nodeIndex     = m_nodeIndex;
}

/**
 * Get the surface node average data source.
 *
 * @param chartableFileName
 *    Name of the chartable file.
 * @param nodeIndices
 *    Indices of the surface nodes.
 */
void
ChartDataSource::getSurfaceNodeAverage(AString& surfaceStructureName,
                                       int32_t& surfaceNumberOfNodes,
                                       std::vector<int32_t>& nodeIndices) const
{
    surfaceStructureName = m_surfaceStructureName;
    surfaceNumberOfNodes = m_surfaceNumberOfNodes;    
    nodeIndices   = m_nodeIndicesAverage;
}

/**
 * Get the volume voxel data source.
 *
 * @param chartableFileName
 *    Name of the chartable file.
 * @param ijk
 *    Indices of the voxel.
 */
void
ChartDataSource::getVolumeVoxel(int64_t ijk[3]) const
{
    ijk[0] = m_voxelIJK[0];
    ijk[1] = m_voxelIJK[1];
    ijk[2] = m_voxelIJK[2];
}

/**
 * Get a description of this object's content.
 * @return String describing this object's content.
 */
AString
ChartDataSource::toString() const
{
    return "ChartDataSource";
}

/**
 * Save information specific to this type of model to the scene.
 *
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    saving the scene.
 *
 * @param instanceName
 *    Name of instance in the scene.
 */
SceneClass*
ChartDataSource::saveToScene(const SceneAttributes* sceneAttributes,
                                          const AString& instanceName)
{
    SceneClass* sceneClass = new SceneClass(instanceName,
                                            "ChartDataSource",
                                            1);
    m_sceneAssistant->saveMembers(sceneAttributes,
                                  sceneClass);
    sceneClass->addPathName("m_chartableFileName",
                            m_chartableFileName);
    
    switch (m_dataSourceMode) {
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_INVALID:
            break;
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_SURFACE_NODE_INDEX:
            break;
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_SURFACE_NODE_INDICES_AVERAGE:
        {
            const int32_t numNodes = static_cast<int32_t>(m_nodeIndicesAverage.size());
            if (numNodes > 0) {
                sceneClass->addIntegerArray("m_nodeIndicesAverage",
                                            &m_nodeIndicesAverage[0],
                                            numNodes);
            }
        }
            break;
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_VOXEL_IJK:
            break;
    }
    
    return sceneClass;
}

/**
 * Restore information specific to the type of model from the scene.
 *
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    restoring the scene.
 *
 * @param sceneClass
 *     sceneClass from which model specific information is obtained.
 */
void
ChartDataSource::restoreFromScene(const SceneAttributes* sceneAttributes,
                                               const SceneClass* sceneClass)
{
    if (sceneClass == NULL) {
        return;
    }
    
    m_sceneAssistant->restoreMembers(sceneAttributes,
                                     sceneClass);
    m_chartableFileName = sceneClass->getPathNameValue("m_chartableFileName");
    
    switch (m_dataSourceMode) {
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_INVALID:
            break;
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_SURFACE_NODE_INDEX:
            break;
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_SURFACE_NODE_INDICES_AVERAGE:
        {
            const SceneClassArray* nodeArray = sceneClass->getClassArray("m_nodeIndicesAverage");
            if (nodeArray != NULL) {
                const int32_t numNodes = nodeArray->getNumberOfArrayElements();
                if (numNodes > 0) {
                    m_nodeIndicesAverage.resize(numNodes);
                    sceneClass->getIntegerArrayValue("m_nodeIndicesAverage",
                                                     &m_nodeIndicesAverage[0],
                                                     numNodes);
                }
            }
        }
            break;
        case ChartDataSourceModeEnum::CHART_DATA_SOURCE_MODE_VOXEL_IJK:
            break;
    }
    
}


