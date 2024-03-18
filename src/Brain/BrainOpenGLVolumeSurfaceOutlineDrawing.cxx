
/*LICENSE_START*/
/*
 *  Copyright (C) 2022 Washington University School of Medicine
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

#define __BRAIN_OPEN_G_L_VOLUME_SURFACE_OUTLINE_DRAWING_DECLARE__
#include "BrainOpenGLVolumeSurfaceOutlineDrawing.h"
#undef __BRAIN_OPEN_G_L_VOLUME_SURFACE_OUTLINE_DRAWING_DECLARE__

#include "CaretAssert.h"
#include "CaretLogger.h"
#include "BrainOpenGLFixedPipeline.h"
#include "ElapsedTimer.h"
#include "GraphicsEngineDataOpenGL.h"
#include "GraphicsPrimitive.h"
#include "GraphicsUtilitiesOpenGL.h"
#include "HistologySlice.h"
#include "MathFunctions.h"
#include "Plane.h"
#include "Surface.h"
#include "SurfaceNodeColoring.h"
#include "SurfacePlaneIntersectionToContour.h"
#include "VolumeMappableInterface.h"
#include "VolumeSurfaceOutlineModel.h"
#include "VolumeSurfaceOutlineModelCacheKey.h"
#include "VolumeSurfaceOutlineSetModel.h"

using namespace caret;


    
/**
 * \class caret::BrainOpenGLVolumeSurfaceOutlineDrawing 
 * \brief Draws volume surface outlines
 * \ingroup Brain
 */

/**
 * Constructor.
 */
BrainOpenGLVolumeSurfaceOutlineDrawing::BrainOpenGLVolumeSurfaceOutlineDrawing()
: CaretObject()
{
    
}

/**
 * Destructor.
 */
BrainOpenGLVolumeSurfaceOutlineDrawing::~BrainOpenGLVolumeSurfaceOutlineDrawing()
{
}

/**
 * Get a description of this object's content.
 * @return String describing this object's content.
 */
AString 
BrainOpenGLVolumeSurfaceOutlineDrawing::toString() const
{
    return "BrainOpenGLVolumeSurfaceOutlineDrawing";
}

/**
 * Draw surface outlines on the volume slices
 *
 * @param histologySlice
 *    The histology slice
 * @param outlineSet
 *    The surface outline set.
 * @param fixedPipelineDrawing
 *    The fixed pipeline drawing.
 * @param useNegativePolygonOffsetFlag
 *    If true, use a negative offset for polygon offset
 */
void
BrainOpenGLVolumeSurfaceOutlineDrawing::drawSurfaceOutline(const HistologySlice* histologySlice,
                                                           VolumeSurfaceOutlineSetModel* outlineSet,
                                                           BrainOpenGLFixedPipeline* fixedPipelineDrawing,
                                                           const bool useNegativePolygonOffsetFlag)
{
    const VolumeMappableInterface* underlayVolume(NULL);
    const Plane stereotaxicPlane(histologySlice->getStereotaxicPlane());
    VolumeSurfaceOutlineModelCacheKey outlineCacheKey(histologySlice,
                                                      histologySlice->getSliceIndex());
    drawSurfaceOutlineCached(histologySlice,
                             underlayVolume,
                             ModelTypeEnum::MODEL_TYPE_HISTOLOGY,
                             stereotaxicPlane,
                             outlineCacheKey,
                             outlineSet,
                             fixedPipelineDrawing,
                             useNegativePolygonOffsetFlag);
}

/**
 * Draw surface outlines on the volume slices
 *
 * @param underlayVolume
 *    The underlay volume
 * @param modelType
 *    Type of model being drawn.
 * @param sliceProjectionType
 *    Projection type (oblique/orthogonal)
 * @param sliceViewPlane
 *    Slice view plane (axial, coronal, parasagittal)
 * @param sliceXYZ
 *    Coordinates of slices
 * @param plane
 *    Plane of the volume slice on which surface outlines are drawn.
 * @param outlineSet
 *    The surface outline set.
 * @param fixedPipelineDrawing
 *    The fixed pipeline drawing.
 * @param useNegativePolygonOffsetFlag
 *    If true, use a negative offset for polygon offset
 */
void
BrainOpenGLVolumeSurfaceOutlineDrawing::drawSurfaceOutline(const VolumeMappableInterface* underlayVolume,
                                                           const ModelTypeEnum::Enum modelType,
                                                           const VolumeSliceProjectionTypeEnum::Enum sliceProjectionType,
                                                           const VolumeSliceViewPlaneEnum::Enum sliceViewPlane,
                                                           const float sliceXYZ[3],
                                                           const Plane& plane,
                                                           const Matrix4x4& displayTransformMatrix,
                                                           const bool displayTransformMatrixValidFlag,
                                                           VolumeSurfaceOutlineSetModel* outlineSet,
                                                           BrainOpenGLFixedPipeline* fixedPipelineDrawing,
                                                           const bool useNegativePolygonOffsetFlag)
{
    bool drawCachedFlag(true);
    
    switch (sliceProjectionType) {
        case VolumeSliceProjectionTypeEnum::VOLUME_SLICE_PROJECTION_ORTHOGONAL:
            break;
        case VolumeSliceProjectionTypeEnum::VOLUME_SLICE_PROJECTION_OBLIQUE:
            break;
        case VolumeSliceProjectionTypeEnum::VOLUME_SLICE_PROJECTION_MPR:
            drawCachedFlag = false;
            break;
        case VolumeSliceProjectionTypeEnum::VOLUME_SLICE_PROJECTION_MPR_THREE:
            drawCachedFlag = false;
            break;
    }
    /*
     * Code still here to allow comparison with
     * previous algorithm
     */
    if (drawCachedFlag) {
        drawSurfaceOutlineCachedOnVolume(underlayVolume,
                                 modelType,
                                 sliceProjectionType,
                                 sliceViewPlane,
                                 sliceXYZ,
                                 plane,
                                 outlineSet,
                                 fixedPipelineDrawing,
                                 useNegativePolygonOffsetFlag);
    }
    else {
        drawSurfaceOutlineNotCached(underlayVolume,
                                    modelType,
                                    plane,
                                    displayTransformMatrix,
                                    displayTransformMatrixValidFlag,
                                    outlineSet,
                                    fixedPipelineDrawing,
                                    useNegativePolygonOffsetFlag);
    }
}

/**
 * Draw surface outlines on the volume slices
 *
 * @param underlayVolume
 *    The underlay volume
 * @param modelType
 *    Type of model being drawn.
 * @param sliceProjectionType
 Type of slice projection
 * @param sliceProjectionType
 *    Type of slice projection
 * @param sliceViewPlane
 *    Slice view plane (axial, coronal, parasagittal)
 * @param sliceXYZ
 *    Coordinates of slices
 * @param plane
 *    Plane of the volume slice on which surface outlines are drawn.
 * @param outlineSet
 *    The surface outline set.
 * @param fixedPipelineDrawing
 *    The fixed pipeline drawing.
 * @param useNegativePolygonOffsetFlag
 *    If true, use a negative offset for polygon offset
 */
void
BrainOpenGLVolumeSurfaceOutlineDrawing::drawSurfaceOutlineCachedOnVolume(const VolumeMappableInterface* underlayVolume,
                                                        const ModelTypeEnum::Enum modelType,
                                                        const VolumeSliceProjectionTypeEnum::Enum sliceProjectionType,
                                                        const VolumeSliceViewPlaneEnum::Enum sliceViewPlane,
                                                        const float sliceXYZ[3],
                                                        const Plane& plane,
                                                        VolumeSurfaceOutlineSetModel* outlineSet,
                                                        BrainOpenGLFixedPipeline* fixedPipelineDrawing,
                                                        const bool useNegativePolygonOffsetFlag)
{
    
    float sliceCoordinate(0.0);
    switch (sliceViewPlane) {
        case VolumeSliceViewPlaneEnum::ALL:
            break;
        case VolumeSliceViewPlaneEnum::AXIAL:
            sliceCoordinate = sliceXYZ[2];
            break;
        case VolumeSliceViewPlaneEnum::CORONAL:
            sliceCoordinate = sliceXYZ[1];
            break;
        case VolumeSliceViewPlaneEnum::PARASAGITTAL:
            sliceCoordinate = sliceXYZ[0];
            break;
    }
    
    /*
     * Key for outline cache
     */
    VolumeSurfaceOutlineModelCacheKey outlineCacheKey(underlayVolume,
                                                      sliceViewPlane,
                                                      sliceCoordinate);
    switch (sliceProjectionType) {
        case VolumeSliceProjectionTypeEnum::VOLUME_SLICE_PROJECTION_OBLIQUE:
            outlineCacheKey = VolumeSurfaceOutlineModelCacheKey(underlayVolume,
                                                                plane);
            break;
        case VolumeSliceProjectionTypeEnum::VOLUME_SLICE_PROJECTION_ORTHOGONAL:
            break;
        case VolumeSliceProjectionTypeEnum::VOLUME_SLICE_PROJECTION_MPR:
            break;
        case VolumeSliceProjectionTypeEnum::VOLUME_SLICE_PROJECTION_MPR_THREE:
            break;
    }
    
    const HistologySlice* histologySlice(NULL);
    drawSurfaceOutlineCached(histologySlice,
                             underlayVolume,
                             modelType,
                             plane,
                             outlineCacheKey,
                             outlineSet,
                             fixedPipelineDrawing,
                             useNegativePolygonOffsetFlag);
}

/**
 * Draw surface outlines on the volume slices
 *
 * @param underlayVolume
 *    The underlay volume
 * @param modelType
 *    Type of model on which outlines are drawn
 * @param plane
 *    Plane of the volume slice on which surface outlines are drawn.
 * @param outlineCacheKey
 *    Key for outline cache
 * @param outlineSet
 *    The surface outline set.
 * @param fixedPipelineDrawing
 *    The fixed pipeline drawing.
 * @param useNegativePolygonOffsetFlag
 *    If true, use a negative offset for polygon offset
 */
void
BrainOpenGLVolumeSurfaceOutlineDrawing::drawSurfaceOutlineCached(const HistologySlice* histologySlice,
                                                                 const VolumeMappableInterface* underlayVolume,
                                                                 const ModelTypeEnum::Enum modelType,
                                                                 const Plane& plane,
                                                                 VolumeSurfaceOutlineModelCacheKey& outlineCacheKey,
                                                                 VolumeSurfaceOutlineSetModel* outlineSet,
                                                                 BrainOpenGLFixedPipeline* fixedPipelineDrawing,
                                                                 const bool useNegativePolygonOffsetFlag)
{
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    
    switch (modelType) {
        case ModelTypeEnum::MODEL_TYPE_CHART:
            break;
        case ModelTypeEnum::MODEL_TYPE_CHART_TWO:
            break;
        case ModelTypeEnum::MODEL_TYPE_INVALID:
            break;
        case ModelTypeEnum::MODEL_TYPE_HISTOLOGY:
            break;
        case  ModelTypeEnum::MODEL_TYPE_MULTI_MEDIA:
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE:
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE_MONTAGE:
            break;
        case ModelTypeEnum::MODEL_TYPE_VOLUME_SLICES:
            break;
        case ModelTypeEnum::MODEL_TYPE_WHOLE_BRAIN:
            /*
             * Enable depth so outlines in front or in back
             * of the slices.  Without this the volume surface
             * outlines "behind" the slices are visible and
             * it looks weird
             */
            glEnable(GL_DEPTH_TEST);
            break;
    }
    
    /*
     * Process each surface outline
     * As of 24 May, "zero" is on top so draw in reverse order
     */
    const int32_t numberOfOutlines = outlineSet->getNumberOfDislayedVolumeSurfaceOutlines();
    for (int32_t io = (numberOfOutlines - 1);
         io >= 0;
         io--) {
        std::vector<GraphicsPrimitive*> contourPrimitives;
        
        VolumeSurfaceOutlineModel* outline = outlineSet->getVolumeSurfaceOutlineModel(io);
        if (outline->isDisplayed()) {
            Surface* surface = outline->getSurface();
            if (surface != NULL) {
                float thicknessPercentage = outline->getThicknessPercentageViewportHeight();
                const float thicknessPixels = outline->getThicknessPixelsObsolete();
                
                /*
                 * Thickness was changed from pixels to percentage viewport height on Feb 02, 2018
                 * If thickness percentage is negative, it was not present in an old
                 * scene so convert pixels to percentage using viewports dimensions
                 */
                if (thicknessPercentage < 0.0f) {
                    thicknessPercentage = GraphicsUtilitiesOpenGL::convertPixelsToPercentageOfViewportHeight(thicknessPixels);
                    if (thicknessPercentage > 0.0f) {
                        outline->setThicknessPercentageViewportHeight(thicknessPercentage);
                    }
                }
                
                if (outline->getOutlineCachePrimitives(histologySlice,
                                                       underlayVolume,
                                                       outlineCacheKey,
                                                       contourPrimitives)) {
                    /* OK, have cached primitives to draw */
                }
                else {
                    CaretColorEnum::Enum outlineColor = CaretColorEnum::BLACK;
                    int32_t colorSourceBrowserTabIndex = -1;
                    
                    VolumeSurfaceOutlineColorOrTabModel* colorOrTabModel = outline->getColorOrTabModel();
                    VolumeSurfaceOutlineColorOrTabModel::Item* selectedColorOrTabItem = colorOrTabModel->getSelectedItem();
                    switch (selectedColorOrTabItem->getItemType()) {
                        case VolumeSurfaceOutlineColorOrTabModel::Item::ITEM_TYPE_BROWSER_TAB:
                            colorSourceBrowserTabIndex = selectedColorOrTabItem->getBrowserTabIndex();
                            outlineColor = CaretColorEnum::CUSTOM;
                            break;
                        case VolumeSurfaceOutlineColorOrTabModel::Item::ITEM_TYPE_COLOR:
                            outlineColor = selectedColorOrTabItem->getColor();
                            break;
                    }
                    const bool surfaceColorFlag = (colorSourceBrowserTabIndex >= 0);
                    
                    float* nodeColoringRGBA = NULL;
                    if (surfaceColorFlag) {
                        nodeColoringRGBA = fixedPipelineDrawing->surfaceNodeColoring->colorSurfaceNodes(NULL,
                                                                                                        surface,
                                                                                                        colorSourceBrowserTabIndex);
                    }
                    
                    const float slicePlaneDepth(outline->getSlicePlaneDepth());
                    createContours(underlayVolume,
                                   surface,
                                   plane,
                                   outlineColor,
                                   nodeColoringRGBA,
                                   thicknessPercentage,
                                   slicePlaneDepth,
                                   contourPrimitives);

                    if (histologySlice != NULL) {
                        projectContoursToHistologySlice(histologySlice,
                                                        contourPrimitives);
                    }
                    
                    outline->setOutlineCachePrimitives(histologySlice,
                                                       underlayVolume,
                                                       outlineCacheKey,
                                                       contourPrimitives);
                }
            }
        }
        
        /**
         * Draw the contours.
         * Note: The primitives are now cached so DO NOT delete them.
         */
        for (auto primitive : contourPrimitives) {
            if (useNegativePolygonOffsetFlag) {
                glPolygonOffset(-1.0, -1.0);
            }
            else {
                glPolygonOffset(1.0, 1.0);
            }
            glEnable(GL_POLYGON_OFFSET_FILL);
            
            GraphicsEngineDataOpenGL::draw(primitive);
            
            glDisable(GL_POLYGON_OFFSET_FILL);
        }
    }
    
    glPopAttrib();
}

/**
 * Project the contours that are in sterotaxic space to plane coordinates on the histology slice
 * @param histologySlice
 *    The histology slice
 * @param contourPrimitives
 *    Primitives containing the contours
 */
void
BrainOpenGLVolumeSurfaceOutlineDrawing::projectContoursToHistologySlice(const HistologySlice* histologySlice,
                                                                        std::vector<GraphicsPrimitive*>& contourPrimitives)
{
    for (GraphicsPrimitive* primitive : contourPrimitives) {
        CaretAssert(primitive);
        const std::vector<float>& coordinates(primitive->getFloatXYZ());
        const int32_t numCoordinates(primitive->getNumberOfVertices());
        std::vector<float> planeCoordinates;
        planeCoordinates.reserve(numCoordinates * 3);
        for (int32_t i = 0; i < numCoordinates; i++) {
            Vector3D planeXYZ;
            histologySlice->stereotaxicXyzToPlaneXyz(&coordinates[i * 3],
                                                     planeXYZ);
            planeCoordinates.insert(planeCoordinates.end(),
                                    &planeXYZ[0], &planeXYZ[0] + 3);
        }
        
        CaretAssert(planeCoordinates.size() == coordinates.size());
        primitive->replaceFloatXYZ(planeCoordinates);
    }
}

/**
 * Draw surface outlines on the volume slices WITHOUT caching
 *
 * @param underlayVolume
 *    The intersection volume (NULL if not drawing on a volume)
 * @param modelType
 *    Type of model being drawn.
 * @param plane
 *    Plane of the volume slice on which surface outlines are drawn.
 * @param outlineSet
 *    The surface outline set.
 * @param fixedPipelineDrawing
 *    The fixed pipeline drawing.
 * @param useNegativePolygonOffsetFlag
 *    If true, use a negative offset for polygon offset
 */
void
BrainOpenGLVolumeSurfaceOutlineDrawing::drawSurfaceOutlineNotCached(const VolumeMappableInterface* underlayVolume,
                                                                    const ModelTypeEnum::Enum modelType,
                                                                    const Plane& plane,
                                                                    const Matrix4x4& displayTransformMatrix,
                                                                    const bool displayTransformMatrixValidFlag,
                                                                    VolumeSurfaceOutlineSetModel* outlineSet,
                                                                    BrainOpenGLFixedPipeline* fixedPipelineDrawing,
                                                                    const bool useNegativePolygonOffsetFlag)
{
    glPushAttrib(GL_ENABLE_BIT);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LIGHTING);
    
    switch (modelType) {
        case ModelTypeEnum::MODEL_TYPE_CHART:
            break;
        case ModelTypeEnum::MODEL_TYPE_CHART_TWO:
            break;
        case ModelTypeEnum::MODEL_TYPE_INVALID:
            break;
        case ModelTypeEnum::MODEL_TYPE_HISTOLOGY:
            break;
        case  ModelTypeEnum::MODEL_TYPE_MULTI_MEDIA:
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE:
            break;
        case ModelTypeEnum::MODEL_TYPE_SURFACE_MONTAGE:
            break;
        case ModelTypeEnum::MODEL_TYPE_VOLUME_SLICES:
            break;
        case ModelTypeEnum::MODEL_TYPE_WHOLE_BRAIN:
            /*
             * Enable depth so outlines in front or in back
             * of the slices.  Without this the volume surface
             * outlines "behind" the slices are visible and
             * it looks weird
             */
            glEnable(GL_DEPTH_TEST);
            break;
    }
    
    /*
     * Process each surface outline
     * As of 24 May, "zero" is on top so draw in reverse order
     */
    const int32_t numberOfOutlines = outlineSet->getNumberOfDislayedVolumeSurfaceOutlines();
    for (int32_t io = (numberOfOutlines - 1);
         io >= 0;
         io--) {
        std::vector<GraphicsPrimitive*> contourPrimitives;
        
        VolumeSurfaceOutlineModel* outline = outlineSet->getVolumeSurfaceOutlineModel(io);
        if (outline->isDisplayed()) {
            Surface* surface = outline->getSurface();
            if (surface != NULL) {
                float thicknessPercentage = outline->getThicknessPercentageViewportHeight();
                const float thicknessPixels = outline->getThicknessPixelsObsolete();
                
                /*
                 * Thickness was changed from pixels to percentage viewport height on Feb 02, 2018
                 * If thickness percentage is negative, it was not present in an old
                 * scene so convert pixels to percentage using viewports dimensions
                 */
                if (thicknessPercentage < 0.0f) {
                    thicknessPercentage = GraphicsUtilitiesOpenGL::convertPixelsToPercentageOfViewportHeight(thicknessPixels);
                    if (thicknessPercentage > 0.0f) {
                        outline->setThicknessPercentageViewportHeight(thicknessPercentage);
                    }
                }
                
                CaretColorEnum::Enum outlineColor = CaretColorEnum::BLACK;
                int32_t colorSourceBrowserTabIndex = -1;
                
                VolumeSurfaceOutlineColorOrTabModel* colorOrTabModel = outline->getColorOrTabModel();
                VolumeSurfaceOutlineColorOrTabModel::Item* selectedColorOrTabItem = colorOrTabModel->getSelectedItem();
                switch (selectedColorOrTabItem->getItemType()) {
                    case VolumeSurfaceOutlineColorOrTabModel::Item::ITEM_TYPE_BROWSER_TAB:
                        colorSourceBrowserTabIndex = selectedColorOrTabItem->getBrowserTabIndex();
                        outlineColor = CaretColorEnum::CUSTOM;
                        break;
                    case VolumeSurfaceOutlineColorOrTabModel::Item::ITEM_TYPE_COLOR:
                        outlineColor = selectedColorOrTabItem->getColor();
                        break;
                }
                const bool surfaceColorFlag = (colorSourceBrowserTabIndex >= 0);
                
                float* nodeColoringRGBA = NULL;
                if (surfaceColorFlag) {
                    nodeColoringRGBA = fixedPipelineDrawing->surfaceNodeColoring->colorSurfaceNodes(NULL,
                                                                                                    surface,
                                                                                                    colorSourceBrowserTabIndex);
                }
                
                const float slicePlaneDepth(outline->getSlicePlaneDepth());
                createContours(underlayVolume,
                               surface,
                               plane,
                               outlineColor,
                               nodeColoringRGBA,
                               thicknessPercentage,
                               slicePlaneDepth,
                               contourPrimitives);
            }
        }
        
        /**
         * Draw the contours.
         */
        for (auto primitive : contourPrimitives) {
            CaretAssert(primitive);
            if (useNegativePolygonOffsetFlag) {
                glPolygonOffset(-1.0, -1.0);
            }
            else {
                glPolygonOffset(1.0, 1.0);
            }
            glEnable(GL_POLYGON_OFFSET_FILL);
            
            if (displayTransformMatrixValidFlag) {
                const int32_t numVerts(primitive->getNumberOfVertices());
                for (int32_t i = 0; i < numVerts; i++) {
                    Vector3D xyz;
                    primitive->getVertexFloatXYZ(i, xyz);
                    displayTransformMatrix.multiplyPoint3(xyz);
                    primitive->replaceVertexFloatXYZ(i, xyz);
                }
            }
                
            GraphicsEngineDataOpenGL::draw(primitive);
            delete primitive;
            
            glDisable(GL_POLYGON_OFFSET_FILL);
        }
    }
    
    glPopAttrib();
}

/**
 * Constructor.
 *
 * @param underlayVolume
 *    The intersection volume (NULL if not drawing on a volume)
 * @param surfaceFile
 *     The surface file.
 * @param plane
 *     Plane intersected with the surface.
 * @param caretColor
 *     Solid coloring or, if value is CUSTOM, use the vertex coloring
 * @param vertexColoringRGBA
 *     The per-vertex coloring if 'caretColor' is CUSTOM
 * @param contourThicknessPercentOfViewportHeight
 *     Thickness for the contour as a percentage of viewport height.
 */

void
BrainOpenGLVolumeSurfaceOutlineDrawing::createContours(const VolumeMappableInterface* underlayVolume,
                                                       const SurfaceFile* surface,
                                                       const Plane& plane,
                                                       const CaretColorEnum::Enum outlineColor,
                                                       const float* nodeColoringRGBA,
                                                       const float thicknessPercentage,
                                                       const float slicePlaneDepth,
                                                       std::vector<GraphicsPrimitive*>& contourPrimitives)
{
    const bool timingFlag(false);
    ElapsedTimer timer;
    if (timingFlag) {
        timer.start();
    }
    if (slicePlaneDepth > 0.0) {
        int32_t numSteps(0);
        float depthStart(0.0);
        float depthStepSize(0.0);
        computeDepthNumStepsAndStepSize(underlayVolume,
                                        slicePlaneDepth,
                                        numSteps,
                                        depthStart,
                                        depthStepSize);
        const bool debugFlag(false);
        if (debugFlag) {
            std::cout << "SlicePlaneDepth: " << slicePlaneDepth
            << " Start: " << depthStart
            << " Steps: " << numSteps
            << " Size: " << depthStepSize << std::endl;
        }
        
        for (int32_t i = 0; i < numSteps; i++) {
            const float depthOffset(depthStart +
                                    depthStepSize * static_cast<float>(i));
            Plane intersectionPlane(plane);
            intersectionPlane.shiftPlane(depthOffset);
            const Plane drawOnPlane(plane);
            if (debugFlag) {
                std::cout << i << "   Intersect Plane: " << intersectionPlane.toString() << std::endl;
                std::cout << i << "   Draw On Plane: " << drawOnPlane.toString() << std::endl;
            }
            SurfacePlaneIntersectionToContour contour(surface,
                                                      intersectionPlane,
                                                      drawOnPlane,
                                                      outlineColor,
                                                      nodeColoringRGBA,
                                                      thicknessPercentage);
            AString errorMessage;
            if ( ! contour.createContours(contourPrimitives,
                                          errorMessage)) {
                CaretLogSevere(errorMessage);
            }
        }
    }
    else {
        SurfacePlaneIntersectionToContour contour(surface,
                                                  plane,
                                                  outlineColor,
                                                  nodeColoringRGBA,
                                                  thicknessPercentage);
        AString errorMessage;
        if ( ! contour.createContours(contourPrimitives,
                                      errorMessage)) {
            CaretLogSevere(errorMessage);
        }
    }
    if (timingFlag) {
        std::cout << "Time to compute contours: " << timer.getElapsedTimeMilliseconds() << "ms" << std::endl;
    }
}

/**
 * Compute the number of steps and step size for slice plane depth
 * @param underlayVolume
 *    The intersection volume (NULL if not drawing on a volume)
 * @param slicePlaneDepth
 *    Slice plane depth set by user
 * @param numStepsOut
 *    Number of steps output
 * @param depthStartOut
 *    Starting depth value output
 * @param depthStepSizeOut
 *    Depth step size output
 */
void
BrainOpenGLVolumeSurfaceOutlineDrawing::computeDepthNumStepsAndStepSize(const VolumeMappableInterface* underlayVolume,
                                                                        const float slicePlaneDepth,
                                                                        int32_t& numStepsOut,
                                                                        float& depthStartOut,
                                                                        float& depthStepSizeOut)
{
    /*
     * Default to normal (infinitely thin slice)
     */
    numStepsOut      = 1;
    depthStartOut    = 0.0;
    depthStepSizeOut = 0.0;
    if (slicePlaneDepth <= 0.0) {
        return ;
    }
    
    /*
     * If volume valid, use 1/2 voxel size for step
     */
    depthStepSizeOut = 0.5;
    if (underlayVolume != NULL) {
        depthStepSizeOut = underlayVolume->getMaximumVoxelSpacing() / 2.0;
        if (depthStepSizeOut <= 0.0) {
            depthStepSizeOut = 0.5;
        }
    }
    
    /*
     * Set number of steps
     */
    numStepsOut = static_cast<int32_t>(std::ceil(slicePlaneDepth / depthStepSizeOut));
    if (numStepsOut <= 1) {
        numStepsOut      = 1;
        depthStartOut    = 0.0;
        depthStepSizeOut = 0.0;
        return;
    }
    
    /*
     * Always make number of steps odd so that we
     * do an intersection on the plane of the displayed
     * histology or volume slice
     */
    if (MathFunctions::isEvenNumber(numStepsOut)) {
        ++numStepsOut;
    }
    
    /*
     * Reset step size using number of steps but use
     * number of steps minus one so that the middle step
     * is centered on the intersection plane.
     *
     * For example: If there are three steps, this places
     * an intersections at:
     * % at minus half depth
     * % at the interection plane
     * % at plus half depth
     */
    CaretAssert(MathFunctions::isOddNumber(numStepsOut));
    CaretAssert(numStepsOut >= 3);
    depthStepSizeOut = (slicePlaneDepth / (numStepsOut - 1));
    depthStartOut = -slicePlaneDepth / 2.0;
}
