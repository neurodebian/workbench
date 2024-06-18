
/*LICENSE_START*/
/*
 *  Copyright (C) 2017 Washington University School of Medicine
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

#define __GRAPHICS_PRIMITIVE_DECLARE__
#include "GraphicsPrimitive.h"
#undef __GRAPHICS_PRIMITIVE_DECLARE__

#include <algorithm>
#include <cmath>
#include <numeric>

#include "BoundingBox.h"
#include "CaretAssert.h"
#include "CaretLogger.h"
#include "DescriptiveStatistics.h"
#include "EventManager.h"
#include "GraphicsEngineDataOpenGL.h"
#include "GraphicsPrimitiveV3f.h"
#include "GraphicsPrimitiveV3fC4f.h"
#include "GraphicsPrimitiveV3fC4ub.h"
#include "GraphicsPrimitiveV3fN3f.h"
#include "GraphicsPrimitiveV3fN3fC4f.h"
#include "GraphicsPrimitiveV3fN3fC4ub.h"
#include "GraphicsPrimitiveV3fT2f.h"
#include "GraphicsPrimitiveV3fT3f.h"
#include "MathFunctions.h"
#include "Matrix4x4Interface.h"

using namespace caret;


    
/**
 * \class caret::GraphicsPrimitive 
 * \brief Abstract class for a graphics primitive.
 * \ingroup Graphics
 *
 */

/**
 * Constructs an instance with the given vertex, normal vector, and color type.
 *
 * @param vertexDataType
 *     Data type of the vertices.
 * @param normalVectorDataType
 *     Data type of the normal vectors.
 * @param colorDataType
 *     Data type of the colors.
 * @param vertexColorType
 *     Type of vertex coloring
 * @param textureSettings
 *     Settings for primitives that use textures
 * @param primitiveType
 *     Type of primitive drawn (triangles, lines, etc.)
 */
GraphicsPrimitive::GraphicsPrimitive(const VertexDataType           vertexDataType,
                                     const NormalVectorDataType     normalVectorDataType,
                                     const ColorDataType            colorDataType,
                                     const VertexColorType          vertexColorType,
                                     const GraphicsTextureSettings& textureSettings,
                                     const PrimitiveType            primitiveType)
: CaretObject(),
EventListenerInterface(),
m_vertexDataType(vertexDataType),
m_normalVectorDataType(normalVectorDataType),
m_colorDataType(colorDataType),
m_vertexColorType(vertexColorType),
m_textureSettings(textureSettings),
m_primitiveType(primitiveType),
m_boundingBoxValid(false)
{
    invalidateVertexMeasurements();
    
    switch (m_textureSettings.getDimensionType()) {
        case GraphicsTextureSettings::DimensionType::NONE:
            break;
        case GraphicsTextureSettings::DimensionType::FLOAT_STR_2D:
        case GraphicsTextureSettings::DimensionType::FLOAT_STR_3D:
            switch (m_textureSettings.getMipMappingType()) {
                case GraphicsTextureSettings::MipMappingType::DISABLED:
                    switch (m_textureSettings.getMinificationFilter()) {
                        case GraphicsTextureMinificationFilterEnum::NEAREST:
                        case GraphicsTextureMinificationFilterEnum::LINEAR:
                            break;
                        case GraphicsTextureMinificationFilterEnum::LINEAR_MIPMAP_LINEAR:
                        case GraphicsTextureMinificationFilterEnum::LINEAR_MIPMAP_NEAREST:
                        case GraphicsTextureMinificationFilterEnum::NEAREST_MIPMAP_LINEAR:
                        case GraphicsTextureMinificationFilterEnum::NEAREST_MIPMAP_NEAREST:
                            CaretLogSevere("Mip Mapping is disabled but the minification factor="
                                           + GraphicsTextureMinificationFilterEnum::toName(m_textureSettings.getMinificationFilter())
                                           + " requires mip mapping.  ");
                            break;
                    }
                    break;
                case GraphicsTextureSettings::MipMappingType::ENABLED:
                    break;
            }
            break;
    }
}

/**
 * Destructor.
 */
GraphicsPrimitive::~GraphicsPrimitive()
{
    EventManager::get()->removeAllEventsFromListener(this);    
}

/**
 * Copy constructor for graphics primitive.
 * @param obj
 *    Object that is copied.
 */
GraphicsPrimitive::GraphicsPrimitive(const GraphicsPrimitive& obj)
: CaretObject(obj),
 EventListenerInterface(),
 m_vertexDataType(obj.m_vertexDataType),
 m_normalVectorDataType(obj.m_normalVectorDataType),
 m_colorDataType(obj.m_colorDataType),
 m_vertexColorType(obj.m_vertexColorType),
 m_primitiveType(obj.m_primitiveType),
 m_boundingBoxValid(false)
{
    invalidateVertexMeasurements();
    this->copyHelperGraphicsPrimitive(obj);
}

/**
 * Helps with copying an object of this type.
 * @param obj
 *    Object that is copied.
 */
void 
GraphicsPrimitive::copyHelperGraphicsPrimitive(const GraphicsPrimitive& obj)
{
    m_boundingBoxValid            = false;
    m_xyz                         = obj.m_xyz;
    m_floatNormalVectorXYZ        = obj.m_floatNormalVectorXYZ;
    m_floatRGBA                   = obj.m_floatRGBA;
    m_unsignedByteRGBA            = obj.m_unsignedByteRGBA;
    m_floatTextureSTR             = obj.m_floatTextureSTR;
    m_lineWidthType               = obj.m_lineWidthType;
    m_lineWidthValue              = obj.m_lineWidthValue;
    m_pointSizeType               = obj.m_pointSizeType;
    m_pointDiameterValue          = obj.m_pointDiameterValue;
    m_polygonalLinePrimitiveRestartIndices     = obj.m_polygonalLinePrimitiveRestartIndices;
    m_sphereSizeType              = obj.m_sphereSizeType;
    m_sphereDiameterValue         = obj.m_sphereDiameterValue;
    m_textureSettings             = obj.m_textureSettings;
    m_arrayIndicesSubsetFirstVertexIndex = obj.m_arrayIndicesSubsetFirstVertexIndex;
    m_arrayIndicesSubsetCount     = obj.m_arrayIndicesSubsetCount;
    m_voxelColorUpdate            = obj.m_voxelColorUpdate;
    invalidateVertexMeasurements();


    m_graphicsEngineDataForOpenGL.reset();
}

/**
 * Request the capacity of the primitive for the given number of vertices.
 * Using this method can help reduce reallocations as data is added 
 * to the primitive.  Note: the primitive data is stored in 
 * std::vector's and this method functions identically to
 * std::vector::reserve().
 *
 * @param numberOfVertices
 *     Number of vertices
 */
void
GraphicsPrimitive::reserveForNumberOfVertices(const int32_t numberOfVertices)
{
    switch (m_vertexDataType) {
        case VertexDataType::FLOAT_XYZ:
            m_xyz.reserve(numberOfVertices * 3);
            break;
    }
    
    switch (m_normalVectorDataType) {
        case NormalVectorDataType::FLOAT_XYZ:
            m_floatNormalVectorXYZ.reserve(numberOfVertices * 3);
            break;
        case NormalVectorDataType::NONE:
            break;
    }
    
    switch (m_colorDataType) {
        case ColorDataType::NONE:
            break;
        case ColorDataType::FLOAT_RGBA:
            m_floatRGBA.reserve(numberOfVertices * 4);
            break;
        case ColorDataType::UNSIGNED_BYTE_RGBA:
            m_unsignedByteRGBA.reserve(numberOfVertices * 4);
            break;
    }
    
    switch (m_textureSettings.getDimensionType()) {
        case GraphicsTextureSettings::DimensionType::NONE:
            break;
        case GraphicsTextureSettings::DimensionType::FLOAT_STR_2D:
            m_floatTextureSTR.reserve(numberOfVertices * 3);
            break;
        case GraphicsTextureSettings::DimensionType::FLOAT_STR_3D:
            m_floatTextureSTR.reserve(numberOfVertices * 3);
            break;
    }

    invalidateVertexMeasurements();
}

/**
 * @return The coordinates usage type hint for graphics system.
 */
GraphicsPrimitive::UsageType
GraphicsPrimitive::getUsageTypeCoordinates() const
{
    return m_usageTypeCoordinates;
}
/**
 * @return The normals usage type hint for graphics system.
 */
GraphicsPrimitive::UsageType
GraphicsPrimitive::getUsageTypeNormals() const
{
    return m_usageTypeNormals;
}

/**
 * @return The colors usage type hint for graphics system.
 */
GraphicsPrimitive::UsageType
GraphicsPrimitive::getUsageTypeColors() const
{
    return m_usageTypeColors;
}

/**
 * @return The texture coordinates usage type hint for graphics system.
 */
GraphicsPrimitive::UsageType
GraphicsPrimitive::getUsageTypeTextureCoordinates() const
{
    return m_usageTypeTextureCoordinates;
}

/**
 * Set the usage type hint for all data types for graphics system.
 * This will override any of the individual (vertex, normal, etc)
 * usage types so they should be called after calling this method.
 *
 * @param usageType
 *     New value for usage type.
 */
void
GraphicsPrimitive::setUsageTypeAll(const UsageType usageType)
{
    m_usageTypeCoordinates = usageType;
    m_usageTypeNormals  = usageType;
    m_usageTypeColors   = usageType;
    m_usageTypeTextureCoordinates  = usageType;
}
/**
 * Set the coordinates usage type hint for graphics system.
 *
 * @param usageType
 *     New value for usage type.
 */
void
GraphicsPrimitive::setUsageTypeCoordinates(const UsageType usageType)
{
    m_usageTypeCoordinates = usageType;
}

/**
 * Set the normals usage type hint for graphics system.
 *
 * @param usageType
 *     New value for usage type.
 */
void
GraphicsPrimitive::setUsageTypeNormals(const UsageType usageType)
{
    m_usageTypeNormals = usageType;
}

/**
 * Set the colors usage type hint for graphics system.
 *
 * @param usageType
 *     New value for usage type.
 */
void
GraphicsPrimitive::setUsageTypeColors(const UsageType usageType)
{
    m_usageTypeColors = usageType;
}

/**
 * Set the texture coordinates usage type hint for graphics system.
 *
 * @param usageType
 *     New value for usage type.
 */
void
GraphicsPrimitive::setUsageTypeTextureCoordinates(const UsageType usageType)
{
    m_usageTypeTextureCoordinates = usageType;
}

/**
 * @return Is this graphics primitive valid.
 * Primitive is valid if all of these conditions are met
 * (A) ((count(xyz) * 3) == (count(rgba) * 4))
 * (B) (count(xyz) > 0)
 * (C) ((count(xyz) == count(normals)) OR (count(normals) == 0)
 */
bool
GraphicsPrimitive::isValid() const
{
    switch (m_releaseInstanceDataMode) {
        case ReleaseInstanceDataMode::COMPLETED:
            /*
             * Instance data has been removed
             */
            return true;
            break;
        case ReleaseInstanceDataMode::DISABLED:
            break;
        case ReleaseInstanceDataMode::ENABLED:
            break;
    }
    
    switch (m_vertexDataType) {
        case VertexDataType::FLOAT_XYZ:
            break;
    }
    
    const uint32_t numXYZ = m_xyz.size();
    if (numXYZ > 0) {
        switch (m_normalVectorDataType) {
            case NormalVectorDataType::NONE:
                break;
            case NormalVectorDataType::FLOAT_XYZ:
            {
                const uint32_t numNormalXYZ = m_floatNormalVectorXYZ.size();
                if (numNormalXYZ > 0) {
                    if (numNormalXYZ != numXYZ) {
                        CaretLogWarning("ERROR: GraphicsPrimitive XYZ size not matching surface normal size.");
                        return false;
                    }
                }
            }
                break;
        }
        
        bool haveRgbaFlag = false;
        uint32_t numColorRGBA = 0;
        switch (m_colorDataType) {
            case ColorDataType::NONE:
                haveRgbaFlag = false;
                break;
            case ColorDataType::FLOAT_RGBA:
                numColorRGBA = m_floatRGBA.size();
                haveRgbaFlag = true;
                break;
            case ColorDataType::UNSIGNED_BYTE_RGBA:
                numColorRGBA = m_unsignedByteRGBA.size();
                haveRgbaFlag = true;
                break;
        }
        if (haveRgbaFlag) {
            if ((numXYZ / 3) != (numColorRGBA / 4)) {
                CaretLogWarning("ERROR: GraphicsPrimitive XYZ size does not match RGBA size");
                return false;
            }
        }
        
        switch (m_vertexColorType) {
            case VertexColorType::NONE:
                if (haveRgbaFlag) {
                    CaretLogWarning("ERROR: GraphicsPrimitive VertexColorType is NONE but have RGBA ColorDataType");
                }
                break;
            case VertexColorType::PER_VERTEX_RGBA:
            case VertexColorType::SOLID_RGBA:
                if ( ! haveRgbaFlag) {
                    CaretLogWarning("ERROR: GraphicsPrimitive VertexColorType is RGBA but have NONE ColorDataType");
                }
                break;
        }
        
        bool haveTextureFlag = false;
        uint32_t numTextureSTR = 0;
        switch (m_textureSettings.getDimensionType()) {
            case GraphicsTextureSettings::DimensionType::NONE:
                break;
            case GraphicsTextureSettings::DimensionType::FLOAT_STR_2D:
            case GraphicsTextureSettings::DimensionType::FLOAT_STR_3D:
                numTextureSTR = m_floatTextureSTR.size();
                haveTextureFlag = true;
                break;
        }
        if (haveTextureFlag) {
            if (numXYZ != numTextureSTR) {
                CaretLogWarning("ERROR: GraphicsPrimitive XYZ size does not match Texture STR size");
                return false;
            }
            
            if (m_textureSettings.getImageBytesPointer() == NULL) {
                CaretLogWarning("ERROR: Texture has NULL pointer to image data");
                return false;
            }
            switch (m_textureSettings.getPixelFormatType()) {
                case GraphicsTextureSettings::PixelFormatType::NONE:
                    CaretLogWarning("ERROR: GraphicsPrimitive has texture but NONE for pixel format type");
                    break;
                case GraphicsTextureSettings::PixelFormatType::BGR:
                    break;
                case GraphicsTextureSettings::PixelFormatType::BGRA:
                    break;
                case GraphicsTextureSettings::PixelFormatType::BGRX:
                    break;
                case GraphicsTextureSettings::PixelFormatType::RGB:
                    break;
                case GraphicsTextureSettings::PixelFormatType::RGBA:
                    break;
            }
            
            switch (m_textureSettings.getPixelOrigin()) {
                case GraphicsTextureSettings::PixelOrigin::NONE:
                    CaretLogWarning("ERROR: GraphicsPrimitive has texture but NONE for pixel origin");
                    break;
                case GraphicsTextureSettings::PixelOrigin::BOTTOM_LEFT:
                    break;
                case GraphicsTextureSettings::PixelOrigin::TOP_LEFT:
                    break;
            }
        }
        
        switch (m_primitiveType) {
            case PrimitiveType::OPENGL_LINE_LOOP:
            case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_BEVEL_JOIN:
            case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_MITER_JOIN:
            case PrimitiveType::POLYGONAL_LINE_LOOP_BEVEL_JOIN:
            case PrimitiveType::POLYGONAL_LINE_LOOP_MITER_JOIN:
                if (numXYZ < 3) {
                    CaretLogWarning("Line loop must have at least 3 vertices.");
                }
                break;
            case PrimitiveType::OPENGL_LINE_STRIP:
            case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_BEVEL_JOIN:
            case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_MITER_JOIN:
            case PrimitiveType::POLYGONAL_LINE_STRIP_BEVEL_JOIN:
            case PrimitiveType::POLYGONAL_LINE_STRIP_MITER_JOIN:
                if (numXYZ < 2) {
                    CaretLogWarning("Line strip must have at least 2 vertices.");
                }
                break;
            case PrimitiveType::MODEL_SPACE_POLYGONAL_LINES:
            case PrimitiveType::OPENGL_LINES:
            case PrimitiveType::POLYGONAL_LINES:
                if (numXYZ < 2) {
                    CaretLogWarning("Lines must have at least 2 vertices.");
                }
                else {
                    const uint32_t extraVertices = numXYZ % 2;
                    if (extraVertices > 0) {
                        CaretLogWarning("Extra vertices for drawing lines ignored.");
                    }
                }
                break;
            case PrimitiveType::OPENGL_POINTS:
                break;
            case PrimitiveType::OPENGL_TRIANGLE_FAN:
                if (numXYZ < 3) {
                    CaretLogWarning("Triangle fan must have at least 3 vertices.");
                }
                break;
            case PrimitiveType::OPENGL_TRIANGLE_STRIP:
                if (numXYZ < 3) {
                    CaretLogWarning("Triangle strip must have at least 3 vertices.");
                }
                break;
            case PrimitiveType::OPENGL_TRIANGLES:
                if (numXYZ < 3) {
                    CaretLogWarning("Triangles must have at least 3 vertices.");
                }
                else {
                    const uint32_t extraVertices = numXYZ % 3;
                    if (extraVertices > 0) {
                        CaretLogWarning("Extra vertices for drawing triangles ignored.");
                    }
                }
                break;
            case PrimitiveType::SPHERES:
                if (numXYZ < 1) {
                    CaretLogWarning("Spheres must have at least 1 vertices.");
                }
                break;
        }

        return true;
    }
    
    return false;
}

/**
 * Print the primitive's data, data may be long!
 */
void
GraphicsPrimitive::print() const
{
    std::cout << toStringPrivate(true) << std::endl;
}


/**
 * Get a description of this object's content.
 * @return String describing this object's content.
 */
AString
GraphicsPrimitive::toString() const
{
    return toStringPrivate(false);
}

/**
 * @return Type of primitive as text.
 */
AString
GraphicsPrimitive::getPrimitiveTypeAsText() const
{
    AString s;
    
    switch (m_primitiveType) {
        case PrimitiveType::OPENGL_LINE_LOOP:
            s = "OpenGL Line Loop";
            break;
        case PrimitiveType::OPENGL_LINE_STRIP:
            s = "OpenGL Line Strip";
            break;
        case PrimitiveType::OPENGL_LINES:
            s = "OpenGL Lines";
            break;
        case PrimitiveType::OPENGL_POINTS:
            s = "OpenGL Points";
            break;
        case PrimitiveType::OPENGL_TRIANGLE_FAN:
            s = "OpenGL Triangle Fan";
            break;
        case PrimitiveType::OPENGL_TRIANGLE_STRIP:
            s = "OpenGL Triangle Strip";
            break;
        case PrimitiveType::OPENGL_TRIANGLES:
            s = "OpenGL Triangles";
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_BEVEL_JOIN:
            s = "Model Space Polygonal Line Loop Bevel Join";
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_MITER_JOIN:
            s = "Model Space Polygonal Line Loop Meter Join";
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_BEVEL_JOIN:
            s = "Model Space Polygonal Line Strip Bevel Join";
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_MITER_JOIN:
            s = "Model Space Polygonal Line Strip Miter Join";
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINES:
            s = "Model Space Polygonal Lines";
            break;
        case PrimitiveType::POLYGONAL_LINE_LOOP_BEVEL_JOIN:
            s = "Polygonal Line Loop Bevel Join";
            break;
        case PrimitiveType::POLYGONAL_LINE_LOOP_MITER_JOIN:
            s = "Polygonal Line Loop Miter Join";
            break;
        case PrimitiveType::POLYGONAL_LINE_STRIP_BEVEL_JOIN:
            s = "Polygonal Line Strip Bevel Join";
            break;
        case PrimitiveType::POLYGONAL_LINE_STRIP_MITER_JOIN:
            s = "Polygonal Line Strip Miter Join";
            break;
        case PrimitiveType::POLYGONAL_LINES:
            s = "Polygonal Lines";
            break;
        case PrimitiveType::SPHERES:
            s = "Spheres";
            break;
    }
    
    return s;
}

/**
 * @return Sphere size type as text for the given size type.
 *
 * @param sizeType
 *     The size type.
 */
AString
GraphicsPrimitive::getSphereSizeTypeAsText(const SphereSizeType sizeType) const
{
    AString s;
    
    switch (sizeType) {
        case SphereSizeType::MILLIMETERS:
            s = "Millimeters";
            break;
    }
    return s;
}

/**
 * @return Vertex color type as text for the given vertex color type.
 *
 * @param vertexColorType
 *     The vertex color type.
 */
AString
GraphicsPrimitive::getVertexColorTypeAsText(const VertexColorType vertexColorType) const
{
    AString s;
    
    switch (vertexColorType) {
        case VertexColorType::NONE:
            s = "None";
            break;
        case VertexColorType::PER_VERTEX_RGBA:
            s = "Per Vertex RGBA";
            break;
        case VertexColorType::SOLID_RGBA:
            s = "Solid RGBA";
            break;
    }
    
    return s;
}

/**
 * @return Point size type as text for the given size type.
 *
 * @param sizeType
 *     The size type.
 */
AString
GraphicsPrimitive::getPointSizeTypeAsText(const PointSizeType sizeType) const
{
    AString s;
    
    switch (sizeType) {
        case PointSizeType::MILLIMETERS:
            s = "Millimeters";
            break;
        case PointSizeType::PERCENTAGE_VIEWPORT_HEIGHT:
            s = "Percentage Viewport Height";
            break;
        case PointSizeType::PIXELS:
            s = "Pixels";
            break;
    }
    return s;
}


/**
 * @return Line width type as text for the given line width type.
 *
 * @param lineWidthType
 *     The line width type.
 */
AString
GraphicsPrimitive::getLineWidthTypeAsText(const LineWidthType lineWidthType) const
{
    AString s;
    
    switch (lineWidthType) {
        case LineWidthType::PERCENTAGE_VIEWPORT_HEIGHT:
            s = "Percentage Viewport Height";
            break;
        case LineWidthType::PIXELS:
            s = "Pixels";
            break;
    }
    return s;
}

/**
 * Convert to string.
 *
 * @param includeAllDataFlag
 *     If true, numerical values are included.  Else, only type of data.
 */
AString
GraphicsPrimitive::toStringPrivate(const bool includeAllDataFlag) const
{
    AString s("Primitive: " + getPrimitiveTypeAsText());
    
    const int32_t numVertices = getNumberOfVertices();
    s.appendWithNewLine("Number of Vertices: " + AString::number(numVertices) + "\n");
    
    switch (m_textureSettings.getDimensionType()) {
        case GraphicsTextureSettings::DimensionType::NONE:
            break;
        case GraphicsTextureSettings::DimensionType::FLOAT_STR_2D:
        case GraphicsTextureSettings::DimensionType::FLOAT_STR_3D:
            s.appendWithNewLine("Texture: " + AString::number(m_floatTextureSTR.size()) + " Float Texture 0.0 to 1.0.  ");
            s.appendWithNewLine("   Width: " + AString::number(m_textureSettings.getImageWidth())
                                + " Height: " + AString::number(m_textureSettings.getImageHeight())
                                + " Slices: " + AString::number(m_textureSettings.getImageSlices()));
            s.appendWithNewLine("Texture Settings");
            s.appendWithNewLine(m_textureSettings.toString());
            break;
    }
    
    bool addLineWidthFlag  = false;
    bool addPointSizeFlag  = false;
    bool addSphereSizeFlag = false;
    
    switch (m_primitiveType) {
        case PrimitiveType::OPENGL_LINE_LOOP:
            addLineWidthFlag = true;
            break;
        case PrimitiveType::OPENGL_LINE_STRIP:
            addLineWidthFlag = true;
            break;
        case PrimitiveType::OPENGL_LINES:
            addLineWidthFlag = true;
            break;
        case PrimitiveType::OPENGL_POINTS:
            addPointSizeFlag = true;
            break;
        case PrimitiveType::OPENGL_TRIANGLE_FAN:
            break;
        case PrimitiveType::OPENGL_TRIANGLE_STRIP:
            break;
        case PrimitiveType::OPENGL_TRIANGLES:
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_MITER_JOIN:
            addLineWidthFlag = true;
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_MITER_JOIN:
            addLineWidthFlag = true;
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINES:
            addLineWidthFlag = true;
            break;
        case PrimitiveType::POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case PrimitiveType::POLYGONAL_LINE_LOOP_MITER_JOIN:
            addLineWidthFlag = true;
            break;
        case PrimitiveType::POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case PrimitiveType::POLYGONAL_LINE_STRIP_MITER_JOIN:
            addLineWidthFlag = true;
            break;
        case PrimitiveType::POLYGONAL_LINES:
            addLineWidthFlag = true;
            break;
        case PrimitiveType::SPHERES:
            addSphereSizeFlag = true;
            break;
    }
    
    if ( ! m_polygonalLinePrimitiveRestartIndices.empty()) {
        s.appendWithNewLine("Polygonal Line Primitive Restart Indices Count: " + AString::number(m_polygonalLinePrimitiveRestartIndices.size()));
        if (includeAllDataFlag) {
            std::vector<int32_t> temp(m_polygonalLinePrimitiveRestartIndices.begin(),
                                      m_polygonalLinePrimitiveRestartIndices.end());
            s.appendWithNewLine("    " + AString::fromNumbers(temp, ", "));
        }
    }
    
    if (addLineWidthFlag) {
        s.appendWithNewLine("Line Width Type: "
                            + getLineWidthTypeAsText(m_lineWidthType)
                            + "; Value: "
                            + AString::number(m_lineWidthValue, 'f', 3));
    }
    
    if (addPointSizeFlag) {
        s.appendWithNewLine("Point Diameter Type: "
                            + getPointSizeTypeAsText(m_pointSizeType)
                            + "; Value: "
                            + AString::number(m_pointDiameterValue, 'f', 3));
    }
    
    if (addSphereSizeFlag) {
        s.appendWithNewLine("Sphere Diameter Type: "
                            + getSphereSizeTypeAsText(m_sphereSizeType)
                            + "; Value: "
                            + AString::number(m_sphereDiameterValue, 'f', 3));
    }
    
    s.appendWithNewLine("Vertex Color Type: "
                        + getVertexColorTypeAsText(m_vertexColorType));
    
    s.append("\n");
    if (includeAllDataFlag) {
        for (int32_t i = 0; i < numVertices; i++) {
            s.append(AString("%1: ").arg(i, 5));
            
            switch (m_vertexDataType) {
                case VertexDataType::FLOAT_XYZ:
                {
                    CaretAssertVectorIndex(m_xyz, i*3 + 2);
                    const float* xyz = &m_xyz[i * 3];
                    s.append(AString("%1, %2, %3").arg(xyz[0], 10, 'f', 3).arg(xyz[1], 10, 'f', 3).arg(xyz[2], 10, 'f', 3));
                }
                    break;
            }
            
            switch (m_normalVectorDataType) {
                case NormalVectorDataType::FLOAT_XYZ:
                {
                    CaretAssertVectorIndex(m_floatNormalVectorXYZ, i*3 + 2);
                    const float* xyz = &m_floatNormalVectorXYZ[i * 3];
                    s.append(AString("   N:%1, %2, %3").arg(xyz[0], 7, 'f', 5).arg(xyz[1], 7, 'f', 5).arg(xyz[2], 7, 'f', 5));
                }
                    break;
                case NormalVectorDataType::NONE:
                    break;
            }
            
            switch (m_colorDataType) {
                case ColorDataType::NONE:
                    break;
                case ColorDataType::FLOAT_RGBA:
                {
                    CaretAssertVectorIndex(m_floatRGBA, i*4 + 3);
                    const float* rgba = &m_floatRGBA[i * 4];
                    s.append(AString("   RGBAf: %1, %2, %3, %4").arg(rgba[0], 5, 'f', 3).arg(rgba[1], 5, 'f', 3).arg(rgba[2], 5, 'f', 3).arg(rgba[3], 5, 'f', 3));
                }
                    break;
                case ColorDataType::UNSIGNED_BYTE_RGBA:
                {
                    CaretAssertVectorIndex(m_unsignedByteRGBA, i*4 + 3);
                    const uint8_t* rgba = &m_unsignedByteRGBA[i * 4];
                    s.append(AString("   RGBAb: %1, %2, %3, %4").arg(rgba[0], 3).arg(rgba[1], 3).arg(rgba[2], 3).arg(rgba[3], 3));
                }
            }
            
            switch (m_textureSettings.getDimensionType()) {
                case GraphicsTextureSettings::DimensionType::NONE:
                    break;
                case GraphicsTextureSettings::DimensionType::FLOAT_STR_2D:
                case GraphicsTextureSettings::DimensionType::FLOAT_STR_3D:
                {
                    CaretAssertVectorIndex(m_floatTextureSTR, i*3 + 2);
                    const float* str = &m_floatTextureSTR[i * 3];
                    s.append(AString("   T:%1, %2, %3").arg(str[0], 7, 'f', 5).arg(str[1], 7, 'f', 5).arg(str[2], 7, 'f', 5));
                }
                    break;
            }
            s.append("\n");
        }
    }
    
    return s;
}

/**
 * Add a vertex and its normal vector
 *
 * @param xyz
 *     Coordinate of vertex.
 * @param normalVector
 *     Normal vector for the vertex.
 * @param rgbaFloat
 *     Float RGBA color components.
 * @param rgbaByte
 *     Byte RGBA color components.
 */
void
GraphicsPrimitive::addVertexProtected(const float xyz[3],
                                      const float normalVector[3],
                                      const float rgbaFloat[4],
                                      const uint8_t rgbaByte[4],
                                      const float textureSTR[3])
{
    m_xyz.insert(m_xyz.end(),
                 xyz, xyz + 3);
    
    switch (m_normalVectorDataType) {
        case NormalVectorDataType::FLOAT_XYZ:
            CaretAssert(normalVector);
            m_floatNormalVectorXYZ.insert(m_floatNormalVectorXYZ.end(),
                                          normalVector,
                                          normalVector + 3);
            break;
        case NormalVectorDataType::NONE:
            break;
    }
    
    switch (m_colorDataType) {
        case ColorDataType::NONE:
            break;
        case ColorDataType::FLOAT_RGBA:
            CaretAssert(rgbaFloat);
            m_floatRGBA.insert(m_floatRGBA.end(),
                               rgbaFloat,
                               rgbaFloat + 4);
            break;
        case ColorDataType::UNSIGNED_BYTE_RGBA:
            CaretAssert(rgbaByte);
            m_unsignedByteRGBA.insert(m_unsignedByteRGBA.end(),
                                      rgbaByte,
                                      rgbaByte + 4);
            break;
    }
    
    switch (m_textureSettings.getDimensionType()) {
        case GraphicsTextureSettings::DimensionType::NONE:
            break;
        case GraphicsTextureSettings::DimensionType::FLOAT_STR_2D:
        case GraphicsTextureSettings::DimensionType::FLOAT_STR_3D:
            CaretAssert(textureSTR);
            m_floatTextureSTR.insert(m_floatTextureSTR.end(),
                                     textureSTR,
                                     textureSTR + 3);
            break;
    }
    
    invalidateVertexMeasurements();

    /*
     * The triangle strip primitve vertces are filled after two
     * vertices are added after requesting a restart
     */
    if (m_triangleStripPrimitiveRestartIndex == getNumberOfVertices()) {
        m_triangleStripPrimitiveRestartIndex = -1;
        fillTriangleStripPrimitiveRestartVertices();
    }
}

/**
 * Get the XYZ coordinate from the given vertex.
 *
 * @param vertexIndex 
 *     Index of the vertex.
 * @param xyzOut
 *     Output containing the XYZ coordinat.
 */
void
GraphicsPrimitive::getVertexFloatXYZ(const int32_t vertexIndex,
                                     float xyzOut[3]) const
{
    const int32_t i3 = vertexIndex * 3;
    CaretAssertVectorIndex(m_xyz, i3 + 2);
    xyzOut[0] = m_xyz[i3];
    xyzOut[1] = m_xyz[i3+1];
    xyzOut[2] = m_xyz[i3+2];
}

/**
 * Replace the existing XYZ coordinates with the given
 * XYZ coordinates.  The existing and new coordinates
 * MUST BE the same size.
 *
 * @param xyz
 *     The new XYZ coordinates.
 */
void
GraphicsPrimitive::replaceFloatXYZ(const std::vector<float>& xyz)
{
    switch (m_releaseInstanceDataMode) {
        case ReleaseInstanceDataMode::COMPLETED:
            {
                const QString msg("All XYZ Data in primitive cannot be replaced with different sized data.  "
                                  "Instance data was removed to save memory.  "
                                  "setReleaseInstanceDataMode() should not be called for this primitive.");
                CaretAssertMessage(0, msg);
                CaretLogSevere(msg);
                return;
            }
            break;
        case ReleaseInstanceDataMode::DISABLED:
            break;
        case ReleaseInstanceDataMode::ENABLED:
            break;
    }
    
    if (xyz.size() == m_xyz.size()) {
        m_xyz = xyz;
        
        if (m_graphicsEngineDataForOpenGL != NULL) {
            m_graphicsEngineDataForOpenGL->invalidateCoordinates();
        }
    }
    else {
        const AString msg("Replacement XYZ must be same size as existing xyz");
        CaretLogWarning(msg);
        CaretAssertMessage(0, msg);
    }
    
    invalidateVertexMeasurements();
}

/**
 * Get the Y-components from the XYZ coordinates
 * @param yComponentsOut
 *    Output containing the Y-components (this method will resize to correct number of elements)
 */
void
GraphicsPrimitive::getFloatYComponents(std::vector<float>& yComponentsOut) const
{
    const int32_t numXYZ = getNumberOfVertices();
    yComponentsOut.resize(numXYZ);
    
    for (int32_t i = 0; i < numXYZ; i++) {
        const int32_t i3 = (i * 3);
        CaretAssertVectorIndex(yComponentsOut, i);
        CaretAssertVectorIndex(m_xyz, i3 + 1);
        yComponentsOut[i] = m_xyz[i3 + 1];
    }
}

/**
 * Set the Y-components from the XYZ coordinates
 * @param yComponents
 *    The Y-components
 */
void
GraphicsPrimitive::setFloatYComponents(const std::vector<float>& yComponents)
{
    const int32_t numXYZ = getNumberOfVertices();
    if (numXYZ != static_cast<int32_t>(yComponents.size())) {
        const QString msg("Number of XYZ="
                          + QString::number(numXYZ)
                          + " but number of input Y-components="
                          + QString::number(yComponents.size()));
        CaretLogSevere(msg);
        CaretAssertMessage(0, msg);
        return;
    }
    
    for (int32_t i = 0; i < numXYZ; i++) {
        const int32_t i3 = (i * 3);
        CaretAssertVectorIndex(yComponents, i);
        CaretAssertVectorIndex(m_xyz, i3 + 1);
        m_xyz[i3 + 1] = yComponents[i];
    }
    
    invalidateVertexMeasurements();
    if (m_graphicsEngineDataForOpenGL != NULL) {
        m_graphicsEngineDataForOpenGL->invalidateCoordinates();
    }
}

/**
 * @return Number of vertices in primitive
 */
int32_t
GraphicsPrimitive::getNumberOfVertices() const
{
    const int32_t numVertices(m_xyz.size() / 3);
    return numVertices;
}

/**
 * Replace the XYZ coordinate at the given index
 *
 * @param vertexIndex
 *     Index of the vertex
 * @param xyz
 *     The new XYZ coordinate
 */
void
GraphicsPrimitive::replaceVertexFloatXYZ(const int32_t vertexIndex,
                                         const float xyz[3])
{
    switch (m_releaseInstanceDataMode) {
        case ReleaseInstanceDataMode::COMPLETED:
        {
            const QString msg("A vertex's XYZ Data in primitive cannot be replaced.  "
                              "Instance data was removed to save memory.  "
                              "setReleaseInstanceDataMode() should not be called for this primitive.");
            CaretAssertMessage(0, msg);
            CaretLogSevere(msg);
            return;
        }
            break;
        case ReleaseInstanceDataMode::DISABLED:
            break;
        case ReleaseInstanceDataMode::ENABLED:
            break;
    }

    const int32_t offset = vertexIndex * 3;
    CaretAssertVectorIndex(m_xyz, offset + 2);
    m_xyz[offset]     = xyz[0];
    m_xyz[offset + 1] = xyz[1];
    m_xyz[offset + 2] = xyz[2];

    if (m_graphicsEngineDataForOpenGL != NULL) {
        m_graphicsEngineDataForOpenGL->invalidateCoordinates();
    }
}

/**
 * Replace the STR coordinate at the given index
 *
 * @param vertexIndex
 *     Index of the vertex
 * @param str
 *     The new STR texture coordinate
 */
void
GraphicsPrimitive::replaceVertexTextureSTR(const int32_t vertexIndex,
                                           const float str[3])
{
    switch (m_releaseInstanceDataMode) {
        case ReleaseInstanceDataMode::COMPLETED:
        {
            const QString msg("A vertex's STR Data in primitive cannot be replaced.  "
                              "Instance data was removed to save memory.  "
                              "setReleaseInstanceDataMode() should not be called for this primitive.");
            CaretAssertMessage(0, msg);
            CaretLogSevere(msg);
            return;
        }
            break;
        case ReleaseInstanceDataMode::DISABLED:
            break;
        case ReleaseInstanceDataMode::ENABLED:
            break;
    }
    
    const int32_t offset = vertexIndex * 3;
    CaretAssertVectorIndex(m_floatTextureSTR, offset + 2);
    m_floatTextureSTR[offset]     = str[0];
    m_floatTextureSTR[offset + 1] = str[1];
    m_floatTextureSTR[offset + 2] = str[2];
    
    if (m_graphicsEngineDataForOpenGL != NULL) {
        m_graphicsEngineDataForOpenGL->invalidateTextureCoordinates();
    }
}


/**
 * Replace the XYZ vertices in this primitive with vertices from other primitive
 * and also transform the vertex using the given matrix.
 *
 * @param primitive
 *     Primitive whose vertices are copied.
 * @param matrix
 *     Matrix used for transforming vertices
 */
void
GraphicsPrimitive::replaceAndTransformVertices(const GraphicsPrimitive* primitive,
                                               const Matrix4x4Interface& matrix)
{
    CaretAssert(primitive);
    const int32_t numVertices = std::min(getNumberOfVertices(),
                                         primitive->getNumberOfVertices());
    for (int32_t i = 0; i < numVertices; i++) {
        const std::vector<float>& primitiveXYZ = primitive->getFloatXYZ();
        const int32_t i3(i * 3);
        float xyz[3] { primitiveXYZ[i3], primitiveXYZ[i3 + 1], primitiveXYZ[i3 + 2] };
        matrix.multiplyPoint3(xyz);
        replaceVertexFloatXYZ(i, xyz);
    }
}


/**
 * Transform all vertices using the given matrix.
 * @param matrix
 *    Matrix for transformation of vertices
 */
void
GraphicsPrimitive::transformVerticesFloatXYZ(const Matrix4x4Interface& matrix)
{
    const float numVertices = getNumberOfVertices();
    for (int32_t i = 0; i < numVertices; i++) {
        const int32_t i3(i * 3);
        matrix.multiplyPoint3(&m_xyz[i3]);
    }
}

/**
 * Get the float RGBA coloring for a vertex.
 *
 * @param vertexIndex
 *     Index of the vertex.
 * @param rgbaOut
 *     Output with RGBA.
 */
void
GraphicsPrimitive::getVertexFloatRGBA(const int32_t vertexIndex,
                                      float rgbaOut[4]) const
{
    const int32_t i4 = vertexIndex * 4;
    switch (m_colorDataType) {
        case ColorDataType::NONE:
            CaretAssert(0);
            break;
        case ColorDataType::FLOAT_RGBA:
            CaretAssertVectorIndex(m_floatRGBA, i4 + 3);
            rgbaOut[0] = m_floatRGBA[i4];
            rgbaOut[1] = m_floatRGBA[i4 + 1];
            rgbaOut[2] = m_floatRGBA[i4 + 2];
            rgbaOut[3] = m_floatRGBA[i4 + 3];
            break;
        case ColorDataType::UNSIGNED_BYTE_RGBA:
            CaretAssertMessage(0, "Getting float RGBA from primitive but coloring type is Byte");
            CaretLogWarning("Getting float RGBA from primitive but coloring type is Byte");
            break;
    }
}

/**
 * Replace the float RGBA coloring for a vertex.
 *
 * @param vertexIndex
 *     Index of the vertex.
 * @param rgba
 *     New RGBA for vertex.
 */
void
GraphicsPrimitive::replaceVertexFloatRGBA(const int32_t vertexIndex,
                                          const float rgba[4])
{
    switch (m_releaseInstanceDataMode) {
        case ReleaseInstanceDataMode::COMPLETED:
        {
            const QString msg("A vertex's RGBA Data in primitive cannot be replaced.  "
                              "Instance data was removed to save memory.  "
                              "setReleaseInstanceDataMode() should not be called for this primitive.");
            CaretAssertMessage(0, msg);
            CaretLogSevere(msg);
            return;
        }
            break;
        case ReleaseInstanceDataMode::DISABLED:
            break;
        case ReleaseInstanceDataMode::ENABLED:
            break;
    }

    const int32_t i4 = vertexIndex * 4;
    switch (m_colorDataType) {
        case ColorDataType::NONE:
            CaretAssert(0);
            break;
        case ColorDataType::FLOAT_RGBA:
            CaretAssertVectorIndex(m_floatRGBA, i4 + 3);
            m_floatRGBA[i4]     = rgba[0];
            m_floatRGBA[i4 + 1] = rgba[1];
            m_floatRGBA[i4 + 2] = rgba[2];
            m_floatRGBA[i4 + 3] = rgba[3];
            break;
        case ColorDataType::UNSIGNED_BYTE_RGBA:
            CaretAssertMessage(0, "Replacing float RGBA in primitive but coloring type is Byte");
            CaretLogWarning("Replacing float RGBA in primitive but coloring type is Byte");
            break;
    }

    if (m_graphicsEngineDataForOpenGL != NULL) {
        m_graphicsEngineDataForOpenGL->invalidateColors();
    }
}

/**
 * Get the byte RGBA coloring for a vertex.
 *
 * @param vertexIndex
 *     Index of the vertex.
 * @param rgbaOut
 *     Output with RGBA.
 */
void
GraphicsPrimitive::getVertexByteRGBA(const int32_t vertexIndex,
                                     uint8_t rgbaOut[4]) const
{
    const int32_t i4 = vertexIndex * 4;
    switch (m_colorDataType) {
        case ColorDataType::NONE:
            CaretAssert(0);
            break;
        case ColorDataType::FLOAT_RGBA:
            CaretAssertMessage(0, "Getting Byte RGBA in primitive but coloring type is Float");
            CaretLogWarning("Getting Byte RGBA in primitive but coloring type is Float");
            break;
        case ColorDataType::UNSIGNED_BYTE_RGBA:
            CaretAssertVectorIndex(m_unsignedByteRGBA, i4 + 3);
            rgbaOut[0] = m_unsignedByteRGBA[i4];
            rgbaOut[1] = m_unsignedByteRGBA[i4 + 1];
            rgbaOut[2] = m_unsignedByteRGBA[i4 + 2];
            rgbaOut[3] = m_unsignedByteRGBA[i4 + 3];
            break;
    }
}

/**
 * Replace the float RGBA coloring for a vertex.
 *
 * @param vertexIndex
 *     Index of the vertex.
 * @param rgba
 *     New RGBA for vertex.
 */
void
GraphicsPrimitive::replaceVertexByteRGBA(const int32_t vertexIndex,
                                         const uint8_t rgba[4])
{
    switch (m_releaseInstanceDataMode) {
        case ReleaseInstanceDataMode::COMPLETED:
        {
            const QString msg("A vertex's RGBA Data in primitive cannot be replaced.  "
                              "Instance data was removed to save memory.  "
                              "setReleaseInstanceDataMode() should not be called for this primitive.");
            CaretAssertMessage(0, msg);
            CaretLogSevere(msg);
            return;
        }
            break;
        case ReleaseInstanceDataMode::DISABLED:
            break;
        case ReleaseInstanceDataMode::ENABLED:
            break;
    }

    const int32_t i4 = vertexIndex * 4;
    switch (m_colorDataType) {
        case ColorDataType::NONE:
            CaretAssert(0);
            break;
        case ColorDataType::FLOAT_RGBA:
            CaretAssertMessage(0, "Replacing Byte RGBA in primitive but coloring type is Float");
            CaretLogWarning("Replacing Byte RGBA in primitive but coloring type is Float");
            break;
        case ColorDataType::UNSIGNED_BYTE_RGBA:
            CaretAssertVectorIndex(m_unsignedByteRGBA, i4 + 3);
            m_unsignedByteRGBA[i4]     = rgba[0];
            m_unsignedByteRGBA[i4 + 1] = rgba[1];
            m_unsignedByteRGBA[i4 + 2] = rgba[2];
            m_unsignedByteRGBA[i4 + 3] = rgba[3];
            break;
    }

    if (m_graphicsEngineDataForOpenGL != NULL) {
        m_graphicsEngineDataForOpenGL->invalidateColors();
    }
}

/**
 * Replace all vertex RGBA with the given solid color RGBA.
 *
 * @param rgba
 *     RGBA values for all vertices
 */
void
GraphicsPrimitive::replaceAllVertexSolidByteRGBA(const uint8_t rgba[4])
{
    switch (m_releaseInstanceDataMode) {
        case ReleaseInstanceDataMode::COMPLETED:
        {
            const QString msg("All RGBA Data in primitive cannot be replaced.  "
                              "Instance data was removed to save memory.  "
                              "setReleaseInstanceDataMode() should not be called for this primitive.");
            CaretAssertMessage(0, msg);
            CaretLogSevere(msg);
            return;
        }
            break;
        case ReleaseInstanceDataMode::DISABLED:
            break;
        case ReleaseInstanceDataMode::ENABLED:
            break;
    }

    switch (m_colorDataType) {
        case ColorDataType::NONE:
            CaretAssert(0);
            break;
        case ColorDataType::FLOAT_RGBA:
            CaretAssertMessage(0, "Replacing Byte RGBA in primitive but coloring type is Float");
            CaretLogWarning("Replacing Byte RGBA in primitive but coloring type is Float");
            break;
        case ColorDataType::UNSIGNED_BYTE_RGBA:
        {
            const int32_t numRGBA = static_cast<int32_t>(m_unsignedByteRGBA.size() / 4);
            for (int32_t i = 0; i < numRGBA; i++) {
                const int32_t i4 = i * 4;
                CaretAssertVectorIndex(m_unsignedByteRGBA, i4 + 3);
                m_unsignedByteRGBA[i4]     = rgba[0];
                m_unsignedByteRGBA[i4 + 1] = rgba[1];
                m_unsignedByteRGBA[i4 + 2] = rgba[2];
                m_unsignedByteRGBA[i4 + 3] = rgba[3];
            }
        }
    }
    
    if (m_graphicsEngineDataForOpenGL != NULL) {
        m_graphicsEngineDataForOpenGL->invalidateColors();
    }
}

/**
 * Replace all vertex RGBA with the given solid color RGBA.
 *
 * @param rgba
 *     RGBA values for all vertices
 */
void
GraphicsPrimitive::replaceAllVertexSolidFloatRGBA(const float rgba[4])
{
    switch (m_releaseInstanceDataMode) {
        case ReleaseInstanceDataMode::COMPLETED:
        {
            const QString msg("All RGBA Data in primitive cannot be replaced.  "
                              "Instance data was removed to save memory.  "
                              "setReleaseInstanceDataMode() should not be called for this primitive.");
            CaretAssertMessage(0, msg);
            CaretLogSevere(msg);
            return;
        }
            break;
        case ReleaseInstanceDataMode::DISABLED:
            break;
        case ReleaseInstanceDataMode::ENABLED:
            break;
    }

    switch (m_colorDataType) {
        case ColorDataType::NONE:
            CaretAssert(0);
            break;
        case ColorDataType::FLOAT_RGBA:
        {
            const int32_t numRGBA = static_cast<int32_t>(m_floatRGBA.size() / 4);
            for (int32_t i = 0; i < numRGBA; i++) {
                const int32_t i4 = i * 4;
                CaretAssertVectorIndex(m_floatRGBA, i4 + 3);
                m_floatRGBA[i4]     = rgba[0];
                m_floatRGBA[i4 + 1] = rgba[1];
                m_floatRGBA[i4 + 2] = rgba[2];
                m_floatRGBA[i4 + 3] = rgba[3];
            }
        }
            break;
        case ColorDataType::UNSIGNED_BYTE_RGBA:
            CaretAssertMessage(0, "Replacing Byte RGBA in primitive but coloring type is Float");
            CaretLogWarning("Replacing Byte RGBA in primitive but coloring type is Float");
            break;
    }

    if (m_graphicsEngineDataForOpenGL != NULL) {
        m_graphicsEngineDataForOpenGL->invalidateColors();
    }
}




/**
 * Get a bounds box for the vertex coordinates.
 *
 * @param boundingBoxOut
 *     Output bounding box.
 * @return
 *     True if bounding box is valid, else false.
 */
bool
GraphicsPrimitive::getVertexBounds(BoundingBox& boundingBoxOut) const
{
    if ( ! m_boundingBoxValid) {
        if (m_boundingBox == NULL) {
            m_boundingBox.reset(new BoundingBox());
        }
        m_boundingBox->resetForUpdate();
        
        int32_t firstVertexIndex(0);
        int32_t numberOfVertices(static_cast<int32_t>(m_xyz.size() / 3));
        int32_t subsetFirstVertex(-1);
        int32_t subsetVertexCount(-1);
        if (getDrawArrayIndicesSubset(subsetFirstVertex,
                                      subsetVertexCount)) {
            firstVertexIndex  = subsetFirstVertex;
            numberOfVertices  = subsetVertexCount;
        }
        
        for (int32_t i = 0; i < numberOfVertices; i++) {
            const int32_t i3 = (firstVertexIndex + i) * 3;
            CaretAssertVectorIndex(m_xyz, i3 + 2);
            m_boundingBox->updateExcludeNanInf(&m_xyz[i3]);
        }
        
        m_boundingBoxValid = m_boundingBox->isValid2D();
    }
    
    boundingBoxOut = *m_boundingBox;
    return m_boundingBoxValid;
}

/**
 * Get indices for drawing a subset of the coordinates in the primitive
 * @param firstVertexIndexOut
 *    Output with index of first vertex to draw
 * @param vertexCountOut
 *    Output with number of vertices to draw starting at firstVertexIndexOut
 * @return
 *    True if the output first vertex and count are valid.
 *    False implies that all vertices should be drawn.
 */
bool
GraphicsPrimitive::getDrawArrayIndicesSubset(int32_t& firstVertexIndexOut,
                                             int32_t& vertexCountOut) const
{
    firstVertexIndexOut = m_arrayIndicesSubsetFirstVertexIndex;
    vertexCountOut      = m_arrayIndicesSubsetCount;
    return ((firstVertexIndexOut >= 0)
            && (vertexCountOut > 0));
}

/**
 * Set indices for drawing a subset of the coordinates in the primitive.
 * Set to negative numbers to disable drawing a subset of vertices and instead
 * draw all of the vertices.
 * @param firstVertexIndex
 *    Index of first vertex to draw
 * @param vertexCount
 *    Number of vertices to draw starting at firstVertexIndexOut
 */
void
GraphicsPrimitive::setDrawArrayIndicesSubset(const int32_t firstVertexIndex,
                               const int32_t vertexCount) const
{
    m_arrayIndicesSubsetFirstVertexIndex = firstVertexIndex;
    m_arrayIndicesSubsetCount       = vertexCount;
    m_boundingBoxValid              = false;
}


/**
 * There may be instances where a primitive should stop and restart at
 * non-consecutive vertices (such as a gap in connected line segment).
 * Typically, this requires creating separate primitives around
 * the gap.
 * 
 * When this method is called, it indicates that any previously
 * added points will not connect to any points after this method
 * is called.  
 *
 * This functionality is similar to that of glPrimitiveRestartIndex()
 * that is not available in versions of OpenGL prior to 3.1.
 *
 * At this time, this functionality if available only for 
 * line strip and line loop.
 */
void
GraphicsPrimitive::addPrimitiveRestart()
{
    bool polygonalLineFlag = false;
    bool triangleStripFlag = false;
    
    switch (m_primitiveType) {
        case PrimitiveType::OPENGL_LINE_LOOP:
            break;
        case PrimitiveType::OPENGL_LINE_STRIP:
            break;
        case PrimitiveType::OPENGL_LINES:
            break;
        case PrimitiveType::OPENGL_POINTS:
            break;
        case PrimitiveType::OPENGL_TRIANGLE_FAN:
            break;
        case PrimitiveType::OPENGL_TRIANGLE_STRIP:
            triangleStripFlag = true;
            break;
        case PrimitiveType::OPENGL_TRIANGLES:
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_MITER_JOIN:
            polygonalLineFlag = true;
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_MITER_JOIN:
            polygonalLineFlag = true;
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINES:
            break;
        case PrimitiveType::POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case PrimitiveType::POLYGONAL_LINE_LOOP_MITER_JOIN:
            polygonalLineFlag = true;
            break;
        case PrimitiveType::POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case PrimitiveType::POLYGONAL_LINE_STRIP_MITER_JOIN:
            polygonalLineFlag = true;
            break;
        case PrimitiveType::POLYGONAL_LINES:
            break;
        case PrimitiveType::SPHERES:
            break;
    }
    
    if (polygonalLineFlag) {
        /*
         * A set is used because a vertex could be added more than one time
         */
        m_polygonalLinePrimitiveRestartIndices.insert(getNumberOfVertices());
    }
    else if (triangleStripFlag) {
        const int32_t numVertices = getNumberOfVertices();
        if (numVertices > 3) {
            /*
             * Add 18 dummy vertices that will be used for the
             * 6 degenerate triangles at a later time
             */
            float float4[4] = { 0.0f, 0.0f, 0.0f, 0.0f };
            uint8_t byte4[4] = { 0, 0, 0, 0 };
            for (int32_t i = 0; i < 18; i++) {
                addVertexProtected(float4, float4, float4, byte4, float4);
            }
            
            /*
             * When the number of vertices becomes the index below,
             * the degenerate triangles will be added.
             */
            m_triangleStripPrimitiveRestartIndex = getNumberOfVertices() + 2;
        }
        else {
            CaretLogSevere("Must have at least three vertices before requesting a triangle strip primitive restart.");
        }
    }
    else {
        CaretLogSevere("Requested primitive restart on not yet supported primitive type: "
                       + getPrimitiveTypeAsText());
    }
}

/**
 * Copy a vertex to another vertex.
 *
 * @param copyFromIndex
 *     Index of the 'copied' vertex.
 * @param copyToIndex
 */
void
GraphicsPrimitive::copyVertex(const int32_t copyFromIndex,
                              const int32_t copyToIndex)
{
    const int32_t from3 = copyFromIndex * 3;
    const int32_t to3   = copyToIndex * 3;
    
    for (int32_t i = 0; i < 3; i++) {
        CaretAssertVectorIndex(m_xyz, from3 + i);
        CaretAssertVectorIndex(m_xyz, to3 + i);
        m_xyz[to3 + i] = m_xyz[from3 + i];
        
        switch (m_normalVectorDataType) {
            case NormalVectorDataType::FLOAT_XYZ:
                CaretAssertVectorIndex(m_floatNormalVectorXYZ, from3 + i);
                CaretAssertVectorIndex(m_floatNormalVectorXYZ, to3 + i);
                m_floatNormalVectorXYZ[to3 + i] = m_floatNormalVectorXYZ[from3 + i];
                break;
            case NormalVectorDataType::NONE:
                break;
        }
        
        switch (m_textureSettings.getDimensionType()) {
            case GraphicsTextureSettings::DimensionType::NONE:
                break;
            case GraphicsTextureSettings::DimensionType::FLOAT_STR_2D:
            case GraphicsTextureSettings::DimensionType::FLOAT_STR_3D:
                CaretAssertVectorIndex(m_floatTextureSTR, from3 + i);
                CaretAssertVectorIndex(m_floatTextureSTR, to3 + i);
                m_floatTextureSTR[to3 + i] = m_floatTextureSTR[from3 + i];
                break;
        }
    }
    
    const int32_t from4 = copyFromIndex * 4;
    const int32_t to4   = copyToIndex * 4;
    
    for (int32_t i = 0; i < 4; i++) {
        switch (m_colorDataType) {
            case ColorDataType::NONE:
                break;
            case ColorDataType::FLOAT_RGBA:
                CaretAssertVectorIndex(m_floatRGBA, from4 + i);
                CaretAssertVectorIndex(m_floatRGBA, to4 + i);
                m_floatRGBA[to4 + i] = m_floatRGBA[from4 + i];
                break;
            case ColorDataType::UNSIGNED_BYTE_RGBA:
                CaretAssertVectorIndex(m_unsignedByteRGBA, from4 + i);
                CaretAssertVectorIndex(m_unsignedByteRGBA, to4 + i);
                m_unsignedByteRGBA[to4 + i] = m_unsignedByteRGBA[from4 + i];
                break;
        }
    }
    
    invalidateVertexMeasurements();
}

/**
 * Fill the degenerate vertices for the triangle strip.
 * Since OpenGL's primitive restart (glPrimitiveRestartIndices())
 * function is not available prior to OpenGL 3.1, add
 * degenerate triangles as explained in the following links.
 *
 * @param newStripVertexOneXYZ
 *     First vertex in new strip that will be added after calling this method.
 * @param newStripVertexTwoXYZ
 *     Second vertex in new strip that will be added after calling this method.
 *
 * @see http://www.gamedev.net/page/resources/_/technical/graphics-programming-and-theory/concatenating-triangle-strips-r1871
 * @see http://en.wikipedia.org/wiki/Triangle_strip
 */
void
GraphicsPrimitive::fillTriangleStripPrimitiveRestartVertices()
{
    const int32_t numVertices = getNumberOfVertices();
    
    /* 
     * Vertices "F" and "F" were added just before the call to addPrimitiveRestart()
     * Vertices "G" and "H" were added after the call to addPrimitiveRestart() and
     * this method was called immediately after "H" was added
     * addPrimitiveRestart() added 18 vertices:
     * DEF + EFF + FFG + FGG + GGH + GHG + GGH + GHI
     *
     * From: https://www.gamedev.net/articles/programming/graphics/concatenating-triangle-strips-r1871
     */
    
    const int32_t vertexH = numVertices - 1;
    const int32_t vertexG = numVertices - 2;
    const int32_t vertexF = numVertices - 21;
    const int32_t vertexE = numVertices - 22;
    int32_t copyToIndex = vertexF + 1;
    
    copyVertex(vertexE, copyToIndex++);
    copyVertex(vertexF, copyToIndex++);
    copyVertex(vertexF, copyToIndex++);
    
    copyVertex(vertexF, copyToIndex++);
    copyVertex(vertexF, copyToIndex++);
    copyVertex(vertexG, copyToIndex++);
    
    copyVertex(vertexF, copyToIndex++);
    copyVertex(vertexG, copyToIndex++);
    copyVertex(vertexG, copyToIndex++);
    
    copyVertex(vertexG, copyToIndex++);
    copyVertex(vertexG, copyToIndex++);
    copyVertex(vertexH, copyToIndex++);
    
    copyVertex(vertexG, copyToIndex++);
    copyVertex(vertexH, copyToIndex++);
    copyVertex(vertexG, copyToIndex++);
    
    copyVertex(vertexG, copyToIndex++);
    copyVertex(vertexG, copyToIndex++);
    copyVertex(vertexH, copyToIndex++);
    
    CaretAssert(copyToIndex == vertexG);
}

/**
 * Get the point diameter.
 *
 * @param sizeTypeOut
 *     Type of sizing.
 * @param pointDiameterOut
 *     Diameter of point.
 */
void
GraphicsPrimitive::getPointDiameter(PointSizeType& sizeTypeOut,
                                    float& pointDiameterOut) const
{
    sizeTypeOut  = m_pointSizeType;
    pointDiameterOut = m_pointDiameterValue;
}

/**
 * Set the point diameter.
 *
 * @param sizeType
 *     Type of sizing.
 * @param pointDiameter
 *     Diameter of point.
 */
void
GraphicsPrimitive::setPointDiameter(const PointSizeType sizeType,
                                    const float pointDiameter) const
{
    m_pointSizeType  = sizeType;
    m_pointDiameterValue = pointDiameter;
}

/**
 * Get the line width.
 *
 * @param lineWidthTypeOut
 *     Type of width.
 * @param lineWidthOut
 *     Width of line.
 */
void
GraphicsPrimitive::getLineWidth(LineWidthType& widthTypeOut,
                                float& lineWidthOut) const
{
    widthTypeOut = m_lineWidthType;
    lineWidthOut = m_lineWidthValue;
}

/**
 * Set the line width.
 *
 * @param lineWidthType
 *     Type of width.
 * @param lineWidth
 *     Width of line.
 */
void
GraphicsPrimitive::setLineWidth(const LineWidthType lineWidthType,
                                const float lineWidth) const
{
    m_lineWidthType  = lineWidthType;
    m_lineWidthValue = lineWidth;
}

/**
 * Get the sphere diameter.
 *
 * @param sizeTypeOut
 *     Type of sizing.
 * @param pointDiameterOut
 *     Diameter of sphere.
 */
void
GraphicsPrimitive::getSphereDiameter(SphereSizeType& sizeTypeOut,
                                     float& sphereDiameterOut) const
{
    sizeTypeOut  = m_sphereSizeType;
    sphereDiameterOut = m_sphereDiameterValue;
}

/**
 * Set the sphere diameter.
 *
 * @param sizeType
 *     Type of sizing.
 * @param pointDiameter
 *     Diameter of sphere.
 */
void
GraphicsPrimitive::setSphereDiameter(const SphereSizeType sizeType,
                                     const float sphereDiameter) const
{
    m_sphereSizeType  = sizeType;
    m_sphereDiameterValue = sphereDiameter;
}

/**
 * @return Number of bytes per pixel for the texture pixel format
 */
int32_t
GraphicsPrimitive::getTexturePixelFormatBytesPerPixel() const
{
    int32_t numBytesPerPixel(0);
    switch (m_textureSettings.getPixelFormatType()) {
        case GraphicsTextureSettings::PixelFormatType::NONE:
            numBytesPerPixel = 0;
            break;
        case GraphicsTextureSettings::PixelFormatType::BGR:
            numBytesPerPixel = 3;
            break;
        case GraphicsTextureSettings::PixelFormatType::BGRA:
            numBytesPerPixel = 4;
            break;
        case GraphicsTextureSettings::PixelFormatType::BGRX:
            numBytesPerPixel = 4;
            break;
        case GraphicsTextureSettings::PixelFormatType::RGB:
            numBytesPerPixel = 3;
            break;
        case GraphicsTextureSettings::PixelFormatType::RGBA:
            numBytesPerPixel = 4;
            break;
    }
    
    return numBytesPerPixel;
}

/**
 * Simplify line types by removing every 'skipVertexCount' vertex.
 * When 'skipVertexCount' is 2, every other point in the line is removed.
 * Nothing is done when 'skipVertexCount' is less than 2 or primitive
 * type is not a line.
 * First and last point are alwyas preserved.
 *
 * @param skipVertexCount
 *     Number of vertices to skip.
 */
void
GraphicsPrimitive::simplfyLines(const int32_t skipVertexCount)
{
    if (skipVertexCount < 2) {
        return;
    }
    
    bool lineTypeFlag(false);
    switch (m_primitiveType) {
        case PrimitiveType::OPENGL_LINE_LOOP:
            lineTypeFlag = true;
            break;
        case PrimitiveType::OPENGL_LINE_STRIP:
            lineTypeFlag = true;
            break;
        case PrimitiveType::OPENGL_LINES:
            break;
        case PrimitiveType::OPENGL_POINTS:
            break;
        case PrimitiveType::OPENGL_TRIANGLE_FAN:
            break;
        case PrimitiveType::OPENGL_TRIANGLE_STRIP:
            break;
        case PrimitiveType::OPENGL_TRIANGLES:
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_MITER_JOIN:
            lineTypeFlag = true;
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_MITER_JOIN:
            lineTypeFlag = true;
            break;
        case PrimitiveType::MODEL_SPACE_POLYGONAL_LINES:
            lineTypeFlag = true;
            break;
        case PrimitiveType::POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case PrimitiveType::POLYGONAL_LINE_LOOP_MITER_JOIN:
            lineTypeFlag = true;
            break;
        case PrimitiveType::POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case PrimitiveType::POLYGONAL_LINE_STRIP_MITER_JOIN:
            lineTypeFlag = true;
            break;
        case PrimitiveType::POLYGONAL_LINES:
            lineTypeFlag = true;
            break;
        case PrimitiveType::SPHERES:
            break;
    }
    
    if ( ! lineTypeFlag) {
        return;
    }
    
    std::vector<float> xyz;
    std::vector<float> normals;
    std::vector<float> rgbaFloat;
    std::vector<uint8_t> rgbaByte;
    std::vector<float> textureSTR;
    
    const int32_t numVertices = getNumberOfVertices();
    xyz.reserve(m_xyz.size());
    normals.reserve(m_floatNormalVectorXYZ.size());
    rgbaFloat.reserve(m_floatRGBA.size());
    rgbaByte.reserve(m_unsignedByteRGBA.size());
    textureSTR.reserve(m_floatTextureSTR.size());
    
    const int32_t lastIndex = numVertices - 1;
    for (int32_t i = 0; i < numVertices; i++) {
        bool addFlag(false);
        if (i == 0) {
            addFlag = true;
        }
        else if (i == lastIndex) {
            addFlag = true;
        }
        else {
            const int remainder = i % skipVertexCount;
            if (remainder == 0) {
                addFlag = true;
            }
        }
        
        if (addFlag) {
            std::cout << "Adding vertex: " << i << std::endl;
            
            const int32_t i3 = i * 3;
            const int32_t i4 = i * 4;
            switch (m_vertexDataType) {
                case VertexDataType::FLOAT_XYZ:
                    xyz.push_back(m_xyz[i3]);
                    xyz.push_back(m_xyz[i3+1]);
                    xyz.push_back(m_xyz[i3+2]);
                    break;
            }
            
            switch (m_normalVectorDataType) {
                case NormalVectorDataType::FLOAT_XYZ:
                    normals.push_back(m_floatNormalVectorXYZ[i3]);
                    normals.push_back(m_floatNormalVectorXYZ[i3+1]);
                    normals.push_back(m_floatNormalVectorXYZ[i3+2]);
                    break;
                case NormalVectorDataType::NONE:
                    break;
            }
            
            switch (m_colorDataType) {
                case ColorDataType::FLOAT_RGBA:
                    rgbaFloat.push_back(m_floatRGBA[i4]);
                    rgbaFloat.push_back(m_floatRGBA[i4+1]);
                    rgbaFloat.push_back(m_floatRGBA[i4+2]);
                    rgbaFloat.push_back(m_floatRGBA[i4+3]);
                    break;
                case ColorDataType::UNSIGNED_BYTE_RGBA:
                    rgbaByte.push_back(m_unsignedByteRGBA[i4]);
                    rgbaByte.push_back(m_unsignedByteRGBA[i4+1]);
                    rgbaByte.push_back(m_unsignedByteRGBA[i4+2]);
                    rgbaByte.push_back(m_unsignedByteRGBA[i4+3]);
                    break;
                case ColorDataType::NONE:
                    break;
            }
            
            switch (m_textureSettings.getDimensionType()) {
                case GraphicsTextureSettings::DimensionType::FLOAT_STR_2D:
                case GraphicsTextureSettings::DimensionType::FLOAT_STR_3D:
                    textureSTR.push_back(m_floatTextureSTR[i3]);
                    textureSTR.push_back(m_floatTextureSTR[i3+1]);
                    textureSTR.push_back(m_floatTextureSTR[i3+2]);
                    break;
                case GraphicsTextureSettings::DimensionType::NONE:
                    break;
            }
        }
    }

    m_xyz = std::move(xyz);
    m_floatNormalVectorXYZ = normals;
    m_floatRGBA = rgbaFloat;
    m_unsignedByteRGBA = rgbaByte;
    m_floatTextureSTR = textureSTR;
}

/**
 * Get the mean and standard deviation for the Y-values
 * @param yMeanOut
 *    Output with mean
 * @param yStandardDeviationOut
 *    Output with standard deviation
 */
void
GraphicsPrimitive::getMeanAndStandardDeviationForY(float& yMeanOut,
                                                   float& yStandardDeviationOut) const
{
    /*
     * Standard deviation is always non-negative so use
     * a negative values as invalid or not-computed
     */
    if (m_yStandardDeviation < 0.0) {
        const int32_t numVertices = getNumberOfVertices();
        if (numVertices <= 0) {
            yMeanOut = 0.0;
            yStandardDeviationOut = 1.0;
            return;
        }
        
        std::vector<float> yValues;
        yValues.resize(numVertices);
        for (int32_t i = 0; i < numVertices; i++) {
            yValues[i] = m_xyz[(i * 3) + 1];
        }
        
        DescriptiveStatistics stats;
        stats.update(yValues);
        m_yMean = stats.getMean();
        m_yStandardDeviation = stats.getPopulationStandardDeviation();
    }
    
    yMeanOut              = m_yMean;
    yStandardDeviationOut = m_yStandardDeviation;
}

/**
 * Apply a new mean and/or deviation to the Y-components
 * @param settings
 *   Settings containing how to modify mean and deviation
 * @param haveNanInfFlagOut
 *   Output: True if Not a Number or Infinity found in the data
 */
void
GraphicsPrimitive::applyNewMeanAndDeviationToYComponents(const GraphicsLineMeanDeviationSettings& settings,
                                                         bool& haveNanInfFlagOut)
{
    haveNanInfFlagOut = false;
    const bool anyModsFlag(settings.m_newMeanEnabled
                           || settings.m_newDeviationEnabled
                           || settings.m_multiplyDeviationEnabled
                           || settings.m_absoluteValueEnabled
                           || settings.m_dataOffsetEnabled);
    if ( ! anyModsFlag) {
        return;
    }
    
    std::vector<float> data;
    getFloatYComponents(data);
    
    for (auto& d : data) {
        if ( ! MathFunctions::isNumeric(d)) {
            haveNanInfFlagOut = true;
            break;
        }
    }
    
    if (haveNanInfFlagOut) {
        applyNewMeanAndDeviationToYComponentsWithNaNs(data,
                                                      settings);
    }
    else {
        applyNewMeanAndDeviationToYComponentsNoNaNs(data,
                                                    settings);
    }
    
    setFloatYComponents(data);
}

/**
 * Apply a new mean and/or deviation to the Y-components that have NaNs or INFs
 * @param data
 *    The data for modification
 * @param settings
 *    Settings containing how to modify mean and deviation
 */
void
GraphicsPrimitive::applyNewMeanAndDeviationToYComponentsWithNaNs(std::vector<float>& data,
                                                                 const GraphicsLineMeanDeviationSettings& settings)
{

    /*
     * Sum of data while excluding invalid numbers
     */
    int32_t numData(0);
    double dataSum(0.0);
    for (auto& d : data) {
        if (MathFunctions::isNumeric(d)) {
            if (settings.m_absoluteValueEnabled) {
                if (d < 0.0) {
                    d = -d;
                }
            }
            dataSum += d;
            numData++;
        }
    }
    
    if (numData <= 0) {
        /*
         * All data is NaN of INF
         */
        return;
    }
    
    /*
     * Compute mean
     */
    const double numValuesDouble(numData);
    const double dataMean(dataSum / numValuesDouble);

    /*
     * Subtract mean from data
     */
    for (auto& d : data) {
        d -= dataMean;
    }

    if (settings.m_multiplyDeviationEnabled
        || settings.m_newDeviationEnabled) {
        /*
         * Compute deviation of data.
         * Note: mean has been already been subtracted from data
         */
        double sumSQ(0.0);
        for (auto& d : data) {
            if (MathFunctions::isNumeric(d)) {
                sumSQ += (d * d);
            }
        }
        const double variance(sumSQ / numValuesDouble);
        double dataDeviation((variance > 0.0)
                             ? (std::sqrt(variance))
                             : 0.0);
        
        double deviationRatio(1.0);
        if (settings.m_multiplyDeviationEnabled) {
            deviationRatio = dataDeviation = settings.m_multiplyDeviationValue;
        }
        else if (settings.m_newDeviationEnabled) {
            /*
             * Data = (Data / dataDeviation) * newDevation
             *      => Data * (newDeviation / dataDeviation);
             */
            deviationRatio = ((dataDeviation != 0.0)
                              ? (settings.m_newDeviationValue / dataDeviation)
                              : settings.m_newDeviationValue);
        }
        else {
            CaretAssert(0);
        }

        for (auto& d : data) {
            d *= deviationRatio;
        }
    }

    double totalOffset(0.0);
    if (settings.m_newMeanEnabled) {
        totalOffset += settings.m_newMeanValue;
    }
    else {
        totalOffset += dataMean;
    }
    if (settings.m_dataOffsetEnabled) {
        totalOffset += settings.m_dataOffsetValue;
    }
    
    for (auto& d : data) {
        d += totalOffset;
    }
}

/**
 * Apply a new mean and/or deviation to the Y-components that DO NOT have NaNs or INFs
 * @param data
 *    The data for modification
 * @param settings
 *    Settings containing how to modify mean and deviation
 */
void
GraphicsPrimitive::applyNewMeanAndDeviationToYComponentsNoNaNs(std::vector<float>& data,
                                                               const GraphicsLineMeanDeviationSettings& settings)
{
    const int32_t num = static_cast<int32_t>(data.size());
    if (num < 1) {
        return;
    }
    const double numValuesDouble(num);
    
    /*
     * Compute mean
     */
    double dataSum(0.0);
    for (auto& d : data) {
        if (settings.m_absoluteValueEnabled) {
            if (d < 0.0) {
                d = -d;
            }
        }
        dataSum += d;
    }
    const double dataMean(dataSum / numValuesDouble);

    /*
     * Subtract mean from data
     */
    for (auto& d : data) {
        d -= dataMean;
    }
    
    if (settings.m_multiplyDeviationEnabled
        || settings.m_newDeviationEnabled) {
        /*
         * Compute deviation of data.
         * Note: mean has been already been subtracted from data
         */
        double sumSQ(0.0);
        for (auto& d : data) {
            sumSQ += (d * d);
        }
        const double variance(sumSQ / numValuesDouble);
        double dataDeviation((variance > 0.0)
                             ? (std::sqrt(variance))
                             : 0.0);
        
        double deviationRatio(1.0);
        if (settings.m_multiplyDeviationEnabled) {
            deviationRatio = dataDeviation = settings.m_multiplyDeviationValue;
        }
        else if (settings.m_newDeviationEnabled) {
            /*
             * Data = (Data / dataDeviation) * newDevation
             *      => Data * (newDeviation / dataDeviation);
             */
            deviationRatio = ((dataDeviation != 0.0)
                              ? (settings.m_newDeviationValue / dataDeviation)
                              : settings.m_newDeviationValue);
        }
        else {
            CaretAssert(0);
        }
        
        for (auto& d : data) {
            d *= deviationRatio;
        }
    }
    
    double totalOffset(0.0);
    if (settings.m_newMeanEnabled) {
        totalOffset += settings.m_newMeanValue;
    }
    else {
        totalOffset += dataMean;
    }
    if (settings.m_dataOffsetEnabled) {
        totalOffset += settings.m_dataOffsetValue;
    }

    for (auto& d : data) {
        d += totalOffset;
    }
}

/**
 * @return A description of the mean/deviation operations.
 */
AString
GraphicsPrimitive::getNewMeanDeviationOperationDescriptionInHtml()
{
    AString txt;
    
    txt.appendWithNewLine("<html>");
    txt.appendWithNewLine("Order of data elements transformation");
    txt.appendWithNewLine("<ol>");
    txt.appendWithNewLine("<li> If <i>Absolute Values</i> is checked, convert data elements to absolute values.");
    txt.appendWithNewLine("<li> Compute Mean of Data.");
    txt.appendWithNewLine("<li> Subtract <i>Data Mean</i> from data elements.");
    txt.appendWithNewLine("<li> If <i>Multiply Deviation</i> is checked, multiply Data Deviation by"
                          "<i>Multiply Deviation Value</i>.");
    txt.appendWithNewLine("<li> if <i>New Deviation</i> is checked, multiply data elements by "
                          "<i>(New Deviation</i> / <i>Data Deviation</i>).");
    txt.appendWithNewLine("<li> If <i>New Mean</i> is checked, add <i>New Mean</i> to data elements.  "
                          "Otherwise, add <i>Data Mean</i> to data elements.");
    txt.appendWithNewLine("<li> If <i>Data Offset</i> is checked, add <i>Data Offset Value</i> to Data.");
    txt.appendWithNewLine("</ol>");
    txt.appendWithNewLine("Note: NaNs are ignored in all computations.");
    txt.appendWithNewLine("</html>");
    
    return txt;
}


/**
 * Get the OpenGL graphics engine data in this instance.
 *
 * @return
 *     OpenGL graphics engine data or NULL if not found.
 */
GraphicsEngineDataOpenGL*
GraphicsPrimitive::getGraphicsEngineDataForOpenGL()
{
    return m_graphicsEngineDataForOpenGL.get();
}

/**
 * Set the OpenGL graphics engine data in this instance.
 *
 * @param graphicsEngineForOpenGL
 *     OpenGLgraphics engine for which graphics engine data is desired.  This value may
 *     NULL to remove graphics engine data for the given graphics engine.
 *     This instance will take ownership of this data and handle deletion of it.
 */
void
GraphicsPrimitive::setGraphicsEngineDataForOpenGL(GraphicsEngineDataOpenGL* graphicsEngineDataForOpenGL)
{
    m_graphicsEngineDataForOpenGL.reset(graphicsEngineDataForOpenGL);
}

/**
 * Receive an event.
 *
 * @param event
 *    An event for which this instance is listening.
 */
void
GraphicsPrimitive::receiveEvent(Event* /*event*/)
{
//    if (event->getEventType() == EventTypeEnum::) {
//        <EVENT_CLASS_NAME*> eventName = dynamic_cast<EVENT_CLASS_NAME*>(event);
//        CaretAssert(eventName);
//
//        event->setEventProcessed();
//    }
}

/**
 * @return A new primitive for XYZ with solid color float RGBA.  Caller is responsible
 * for deleting the returned pointer.
 *
 * @param primitiveType
 *     Type of primitive drawn (triangles, lines, etc.)
 */
GraphicsPrimitiveV3f*
GraphicsPrimitive::newPrimitiveV3f(const GraphicsPrimitive::PrimitiveType primitiveType,
                                 const float floatRGBA[4])
{
    GraphicsPrimitiveV3f* primitive = new GraphicsPrimitiveV3f(primitiveType,
                                                               floatRGBA);
    return primitive;
}

/**
 * @return A new primitive for XYZ with solid color float RGBA.  Caller is responsible
 * for deleting the returned pointer.
 *
 * @param primitiveType
 *     Type of primitive drawn (triangles, lines, etc.)
 */
GraphicsPrimitiveV3fN3f*
GraphicsPrimitive::newPrimitiveV3fN3f(const GraphicsPrimitive::PrimitiveType primitiveType,
                                      const uint8_t unsignedByteRGBA[4])
{
    GraphicsPrimitiveV3fN3f* primitive = new GraphicsPrimitiveV3fN3f(primitiveType,
                                                                     unsignedByteRGBA);
    return primitive;
}

/**
 * @return A new primitive for XYZ with normals and RGBS.  Caller is responsible
 * for deleting the returned pointer.
 *
 * @param primitiveType
 *     Type of primitive drawn (triangles, lines, etc.)
 */
GraphicsPrimitiveV3fN3fC4f*
GraphicsPrimitive::newPrimitiveV3fN3fC4f(const GraphicsPrimitive::PrimitiveType primitiveType)
{
    GraphicsPrimitiveV3fN3fC4f* primitive = new GraphicsPrimitiveV3fN3fC4f(primitiveType);
    return primitive;
}

/**
 * @return A new primitive for XYZ with normals and RGBS.  Caller is responsible
 * for deleting the returned pointer.
 *
 * @param primitiveType
 *     Type of primitive drawn (triangles, lines, etc.)
 */
GraphicsPrimitiveV3fN3fC4ub*
GraphicsPrimitive::newPrimitiveV3fN3fC4ub(const GraphicsPrimitive::PrimitiveType primitiveType)
{
    GraphicsPrimitiveV3fN3fC4ub* primitive = new GraphicsPrimitiveV3fN3fC4ub(primitiveType);
    return primitive;
}

/**
 * @return A new primitive for XYZ with solid color unsigned byte RGBA.  Caller is responsible
 * for deleting the returned pointer.
 *
 * @param primitiveType
 *     Type of primitive drawn (triangles, lines, etc.)
 */
GraphicsPrimitiveV3f*
GraphicsPrimitive::newPrimitiveV3f(const GraphicsPrimitive::PrimitiveType primitiveType,
                                 const uint8_t unsignedByteRGBA[4])
{
    GraphicsPrimitiveV3f* primitive = new GraphicsPrimitiveV3f(primitiveType,
                                                               unsignedByteRGBA);
    return primitive;
}

/**
 * @return A new primitive for XYZ with float RGBA.  Caller is responsible
 * for deleting the returned pointer.
 *
 * @param primitiveType
 *     Type of primitive drawn (triangles, lines, etc.)
 */
GraphicsPrimitiveV3fC4f*
GraphicsPrimitive::newPrimitiveV3fC4f(const GraphicsPrimitive::PrimitiveType primitiveType)
{
    GraphicsPrimitiveV3fC4f* primitive = new GraphicsPrimitiveV3fC4f(primitiveType);
    return primitive;
}

/**
 * @return A new primitive for XYZ with unsigned byte RGBA.  Caller is responsible
 * for deleting the returned pointer.
 *
 * @param primitiveType
 *     Type of primitive drawn (triangles, lines, etc.)
 */
GraphicsPrimitiveV3fC4ub*
GraphicsPrimitive::newPrimitiveV3fC4ub(const GraphicsPrimitive::PrimitiveType primitiveType)
{
    GraphicsPrimitiveV3fC4ub* primitive = new GraphicsPrimitiveV3fC4ub(primitiveType);
    return primitive;
}

/**
 * @return A new primitive with XYZ and texture ST.  Caller is responsible for
 * deleting the returned pointer.
 *
 * @param primitiveType
 *     Type of primitive drawn (triangles, lines, etc.)
 * @param textureSettings
 *     Settings for textures
 */
GraphicsPrimitiveV3fT2f*
GraphicsPrimitive::newPrimitiveV3fT2f(const GraphicsPrimitive::PrimitiveType primitiveType,
                                      const GraphicsTextureSettings& textureSettings)
{
    GraphicsPrimitiveV3fT2f* primitive = new GraphicsPrimitiveV3fT2f(primitiveType,
                                                                     textureSettings);
    return primitive;
}

/**
 * @return A new primitive with XYZ and texture STR.  Caller is responsible for
 * deleting the returned pointer.
 *
 * @param primitiveType
 *     Type of primitive drawn (triangles, lines, etc.)
 * @param textureSettings
 *     Settings for textures
 */
GraphicsPrimitiveV3fT3f*
GraphicsPrimitive::newPrimitiveV3fT3f(const GraphicsPrimitive::PrimitiveType primitiveType,
                                      const GraphicsTextureSettings& textureSettings)
{
    GraphicsPrimitiveV3fT3f* primitive = new GraphicsPrimitiveV3fT3f(primitiveType,
                                                                     textureSettings);
    return primitive;
}

/**
 * Called after buffers have been loaded and may release instance data.
 * If NOT called after buffers have been loaded, bad things may happen.
 */
void
GraphicsPrimitive::setOpenGLBuffersHaveBeenLoadedByGraphicsEngine()
{
    switch (m_releaseInstanceDataMode) {
        case ReleaseInstanceDataMode::COMPLETED:
            return;
            break;
        case ReleaseInstanceDataMode::DISABLED:
            return;
            break;
        case ReleaseInstanceDataMode::ENABLED:
        {
            /*
             * Note: calling ".clear()" on a vector will set the size
             * of the vector to zero but does not free the memory.  Instead,
             * use swap() which causes deallocation of the memory.
             * See http://www.cplusplus.com/reference/vector/vector/clear/
             */
            std::vector<float>().swap(m_xyz);
            std::vector<float>().swap(m_floatRGBA);
            std::vector<float>().swap(m_floatNormalVectorXYZ);
            std::vector<uint8_t>().swap(m_unsignedByteRGBA);
            std::vector<float>().swap(m_floatTextureSTR);
            
            m_releaseInstanceDataMode = ReleaseInstanceDataMode::COMPLETED;
        }
            break;
    }
}

/**
 * Invalidate vertex measurements
 */
void
GraphicsPrimitive::invalidateVertexMeasurements()
{
    m_boundingBoxValid  = false;
    m_yMean              = 0.0;
    m_yStandardDeviation = -1.0;
}

/**
 * @return The voxel color update
 * @param mapIndex
 *    The map index
 */
const VoxelColorUpdate*
GraphicsPrimitive::getVoxelColorUpdate() const
{
    return &m_voxelColorUpdate;
}

/**
 * Set the voxel color update
 */
void
GraphicsPrimitive::setVoxelColorUpdate(const VoxelColorUpdate& voxelColorUpdate)
{
    m_voxelColorUpdate = voxelColorUpdate;
}

/**
 * Reset the voxel color update
 */
void
GraphicsPrimitive::resetVoxelColorUpdate()
{
    m_voxelColorUpdate.clear();
}



