
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

#include "BoundingBox.h"
#include "CaretAssert.h"
#include "CaretLogger.h"
#include "EventManager.h"
#include "GraphicsEngineDataOpenGL.h"
#include "GraphicsPrimitiveV3f.h"
#include "GraphicsPrimitiveV3fC4f.h"
#include "GraphicsPrimitiveV3fC4ub.h"
#include "GraphicsPrimitiveV3fN3f.h"
#include "GraphicsPrimitiveV3fT3f.h"

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
 * @param textureDataType
 *     Data type of texture coordinates.
 * @param primitiveType
 *     Type of primitive drawn (triangles, lines, etc.)
 */
GraphicsPrimitive::GraphicsPrimitive(const VertexDataType       vertexDataType,
                                     const NormalVectorDataType normalVectorDataType,
                                     const ColorDataType        colorDataType,
                                     const VertexColorType      vertexColorType,
                                     const TextureDataType      textureDataType,
                                     const PrimitiveType        primitiveType)
: CaretObject(),
 EventListenerInterface(),
 m_vertexDataType(vertexDataType),
 m_normalVectorDataType(normalVectorDataType),
 m_colorDataType(colorDataType),
 m_vertexColorType(vertexColorType),
 m_textureDataType(textureDataType),
 m_primitiveType(primitiveType),
 m_boundingBoxValid(false)
{
    
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
 m_textureDataType(obj.m_textureDataType),
 m_primitiveType(obj.m_primitiveType),
 m_boundingBoxValid(false)
{
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
    m_textureImageBytesRGBA       = obj.m_textureImageBytesRGBA;
    m_textureImageWidth           = obj.m_textureImageWidth;
    m_textureImageHeight          = obj.m_textureImageHeight;

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
    
    switch (m_textureDataType) {
        case GraphicsPrimitive::TextureDataType::NONE:
            break;
        case GraphicsPrimitive::TextureDataType::FLOAT_STR:
            m_floatTextureSTR.reserve(numberOfVertices * 3);
            break;
    }
    
    m_boundingBoxValid = false;
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
        switch (m_textureDataType) {
            case TextureDataType::NONE:
                break;
            case TextureDataType::FLOAT_STR:
                numTextureSTR = m_floatTextureSTR.size();
                haveTextureFlag = true;
                break;
        }
        if (haveTextureFlag) {
            if (numXYZ != numTextureSTR) {
                CaretLogWarning("ERROR: GraphicsPrimitive XYZ size does not match Texture STR size");
                return false;
            }
            
            if (m_textureImageBytesRGBA.empty()) {
                CaretLogWarning("ERROR: GraphicsPrimitive has invalid texture data");
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
    
    switch (m_textureDataType) {
        case TextureDataType::NONE:
            break;
        case TextureDataType::FLOAT_STR:
            s.appendWithNewLine("Texture: " + AString::number(m_floatTextureSTR.size()) + " Float Texture 0.0 to 1.0.  ");
            s.appendWithNewLine("   Width: " + AString::number(m_textureImageWidth)
                                + " Height: " + AString::number(m_textureImageHeight));
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
            
            switch (m_textureDataType) {
                case TextureDataType::NONE:
                    break;
                case TextureDataType::FLOAT_STR:
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
    
    switch (m_textureDataType) {
        case TextureDataType::NONE:
            break;
        case TextureDataType::FLOAT_STR:
            CaretAssert(textureSTR);
            m_floatTextureSTR.insert(m_floatTextureSTR.end(),
                                     textureSTR,
                                     textureSTR + 3);
            break;
    }
    
    m_boundingBoxValid = false;
    
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
    
    m_boundingBoxValid = false;
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
        
        const int32_t numberOfVertices = static_cast<int32_t>(m_xyz.size() / 3);
        for (int32_t i = 0; i < numberOfVertices; i++) {
            const int32_t i3 = i * 3;
            CaretAssertVectorIndex(m_xyz, i3 + 2);
            m_boundingBox->update(&m_xyz[i3]);
        }
        
        m_boundingBoxValid = true;
    }
    
    boundingBoxOut = *m_boundingBox;
    return m_boundingBoxValid;
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
        
        switch (m_textureDataType) {
            case TextureDataType::NONE:
                break;
            case TextureDataType::FLOAT_STR:
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
    
    m_boundingBoxValid = false;
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
 * Set the image for the texture.
 *
 * @param imageBytesRGBA
 *     Bytes containing the image data.  4 bytes per pixel.
 * @param imageWidth
 *     Width of the actual image.
 * @param imageHeight
 *     Height of the image.
 */
void
GraphicsPrimitive::setTextureImage(const uint8_t* imageBytesRGBA,
                                   const int32_t imageWidth,
                                   const int32_t imageHeight)
{
    m_textureImageBytesRGBA.clear();
    m_textureImageWidth  = -1;
    m_textureImageHeight = -1;
    
    const int32_t numBytes = imageWidth * imageHeight * 4;
    if (numBytes > 0) {
        m_textureImageBytesRGBA.reserve(numBytes);
        m_textureImageWidth  = imageWidth;
        m_textureImageHeight = imageHeight;
        m_textureImageBytesRGBA.insert(m_textureImageBytesRGBA.end(),
                                       imageBytesRGBA, imageBytesRGBA + numBytes);
    }
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
            
            switch (m_textureDataType) {
                case TextureDataType::FLOAT_STR:
                    textureSTR.push_back(m_floatTextureSTR[i3]);
                    textureSTR.push_back(m_floatTextureSTR[i3+1]);
                    textureSTR.push_back(m_floatTextureSTR[i3+2]);
                    break;
                case TextureDataType::NONE:
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

GraphicsPrimitiveV3fT3f*
GraphicsPrimitive::newPrimitiveV3fT3f(const GraphicsPrimitive::PrimitiveType primitiveType,
                                                   const uint8_t* imageBytesRGBA,
                                                   const int32_t imageWidth,
                                                   const int32_t imageHeight)
{
    GraphicsPrimitiveV3fT3f* primitive = new GraphicsPrimitiveV3fT3f(primitiveType,
                                                                     imageBytesRGBA,
                                                                     imageWidth,
                                                                     imageHeight);
    return primitive;
}

