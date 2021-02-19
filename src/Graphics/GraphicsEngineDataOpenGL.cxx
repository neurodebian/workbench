
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

#define __GRAPHICS_ENGINE_DATA_OPEN_G_L_DECLARE__
#include "GraphicsEngineDataOpenGL.h"
#undef __GRAPHICS_ENGINE_DATA_OPEN_G_L_DECLARE__

#include "CaretAssert.h"
#include "CaretLogger.h"
#include "DeveloperFlagsEnum.h"
#include "EventGraphicsOpenGLCreateBufferObject.h"
#include "EventGraphicsOpenGLCreateTextureName.h"
#include "EventOpenGLObjectToWindowTransform.h"
#include "EventManager.h"
#include "GraphicsOpenGLBufferObject.h"
#include "GraphicsOpenGLError.h"
#include "GraphicsOpenGLPolylineTriangles.h"
#include "GraphicsOpenGLTextureName.h"
#include "GraphicsPrimitive.h"
#include "GraphicsPrimitiveSelectionHelper.h"
#include "GraphicsPrimitiveV3f.h"
#include "GraphicsShape.h"
#include "GraphicsUtilitiesOpenGL.h"
#include "Matrix4x4.h"

using namespace caret;


    
/**
 * \class caret::GraphicsEngineDataOpenGL 
 * \brief Contains data the OpenGL Graphics Engine may associate with primtive data.
 * \ingroup Graphics
 */

/**
 * Constructor.
 */
GraphicsEngineDataOpenGL::GraphicsEngineDataOpenGL()
: GraphicsEngineData()
{
    
}

/**
 * Destructor.
 */
GraphicsEngineDataOpenGL::~GraphicsEngineDataOpenGL()
{
    deleteBuffers();
}

/**
 * Delete the OpenGL buffer objects.
 */
void
GraphicsEngineDataOpenGL::deleteBuffers()
{
    if (m_textureImageDataName != NULL) {
        delete m_textureImageDataName;
        m_textureImageDataName = NULL;
    }
}

/**
 * Invalidate the coordinates after they have
 * changed in the graphics primitive.
 */
void
GraphicsEngineDataOpenGL::invalidateCoordinates()
{
    m_reloadCoordinatesFlag = true;
}

/**
 * Invalidate the colors after they have
 * changed in the graphics primitive.
 */
void
GraphicsEngineDataOpenGL::invalidateColors()
{
    m_reloadColorsFlag = true;
}

/**
 * Get the OpenGL Buffer Usage Hint from the primitive.
 *
 * @param primitiveUsageType
 *     The graphics primitive usage type
 * @return
 *     The OpenGL Buffer Usage Hint.
 */
GLenum
GraphicsEngineDataOpenGL::getOpenGLBufferUsageHint(const GraphicsPrimitive::UsageType primitiveUsageType) const
{
    GLenum usageHint = GL_STREAM_DRAW;
    switch (primitiveUsageType) {
        case GraphicsPrimitive::UsageType::MODIFIED_ONCE_DRAWN_FEW_TIMES:
            usageHint = GL_STREAM_DRAW;
            break;
        case GraphicsPrimitive::UsageType::MODIFIED_ONCE_DRAWN_MANY_TIMES:
            usageHint = GL_STATIC_DRAW;
            break;
        case GraphicsPrimitive::UsageType::MODIFIED_MANY_DRAWN_MANY_TIMES:
            usageHint = GL_DYNAMIC_DRAW;
            break;
    }
    
    return usageHint;
}

/**
 * Load the coordinate buffer.
 *
 * @param primitive
 *     The graphics primitive that will be drawn.
 */
void
GraphicsEngineDataOpenGL::loadCoordinateBuffer(GraphicsPrimitive* primitive)
{
    CaretAssert(primitive);
    
    GLenum usageHint = getOpenGLBufferUsageHint(primitive->getUsageTypeCoordinates());
    
    
    GLsizei coordinateCount = 0;
    switch (primitive->m_vertexDataType) {
        case GraphicsPrimitive::VertexDataType::FLOAT_XYZ:
            m_coordinateDataType = GL_FLOAT;
            m_coordinatesPerVertex = 3; /* X, Y, Z */
            
            coordinateCount = primitive->m_xyz.size();
            const GLuint xyzSizeBytes = coordinateCount * sizeof(float);
            CaretAssert(xyzSizeBytes > 0);
            const GLvoid* xyzDataPointer = (const GLvoid*)&primitive->m_xyz[0];
            
            /*
             * Coordinate buffer may have been created.
             * Coordinates may change but the number of coordinates does not change.
             */
            if (m_coordinateBufferObject == NULL) {
                EventGraphicsOpenGLCreateBufferObject createEvent;
                EventManager::get()->sendEvent(createEvent.getPointer());
                m_coordinateBufferObject.reset(createEvent.getOpenGLBufferObject());
                CaretAssert(m_coordinateBufferObject);
            }
            
            CaretAssert(m_coordinateBufferObject->getBufferObjectName());
            
            glBindBuffer(GL_ARRAY_BUFFER,
                         m_coordinateBufferObject->getBufferObjectName());
            glBufferData(GL_ARRAY_BUFFER,
                         xyzSizeBytes,
                         xyzDataPointer,
                         usageHint);
            break;
    }
    
    if (m_coordinatesPerVertex > 0) {
        m_arrayIndicesCount = coordinateCount / m_coordinatesPerVertex;
    }
    else {
        m_arrayIndicesCount = 0;
    }
    
    m_reloadCoordinatesFlag = false;
}

/**
 * Load the normal vector buffer.
 *
 * @param primitive
 *     The graphics primitive that will be drawn.
 */
void
GraphicsEngineDataOpenGL::loadNormalVectorBuffer(GraphicsPrimitive* primitive)
{
    CaretAssert(primitive);
    
    GLenum usageHint = getOpenGLBufferUsageHint(primitive->getUsageTypeNormals());
    
    
    switch (primitive->m_normalVectorDataType) {
        case GraphicsPrimitive::NormalVectorDataType::NONE:
            break;
        case GraphicsPrimitive::NormalVectorDataType::FLOAT_XYZ:
        {
            m_normalVectorDataType = GL_FLOAT;
            const GLuint normalSizeBytes = primitive->m_floatNormalVectorXYZ.size() * sizeof(float);
            CaretAssert(normalSizeBytes > 0);
            const GLvoid* normalDataPointer = (const GLvoid*)&primitive->m_floatNormalVectorXYZ[0];
            
            EventGraphicsOpenGLCreateBufferObject createEvent;
            EventManager::get()->sendEvent(createEvent.getPointer());
            m_normalVectorBufferObject.reset(createEvent.getOpenGLBufferObject());
            CaretAssert(m_normalVectorBufferObject->getBufferObjectName());
            
            glBindBuffer(GL_ARRAY_BUFFER,
                         m_normalVectorBufferObject->getBufferObjectName());
            glBufferData(GL_ARRAY_BUFFER,
                         normalSizeBytes,
                         normalDataPointer,
                         usageHint);
        }
            break;
    }
}

/**
 * Load the color buffer.
 * @param primitive
 *     The graphics primitive that will be drawn.
 */
void
GraphicsEngineDataOpenGL::loadColorBuffer(GraphicsPrimitive* primitive)
{
    CaretAssert(primitive);
    
    GLenum usageHint = getOpenGLBufferUsageHint(primitive->getUsageTypeColors());
    
    switch (primitive->m_colorDataType) {
        case GraphicsPrimitive::ColorDataType::NONE:
            break;
        case GraphicsPrimitive::ColorDataType::FLOAT_RGBA:
        {
            EventGraphicsOpenGLCreateBufferObject createEvent;
            EventManager::get()->sendEvent(createEvent.getPointer());
            m_colorBufferObject.reset(createEvent.getOpenGLBufferObject());
            CaretAssert(m_colorBufferObject->getBufferObjectName());
            
            m_componentsPerColor = 4;
            m_colorDataType = GL_FLOAT;
            
            const GLuint colorSizeBytes = primitive->m_floatRGBA.size() * sizeof(float);
            CaretAssert(colorSizeBytes > 0);
            const GLvoid* colorDataPointer = (const GLvoid*)&primitive->m_floatRGBA[0];
            
            glBindBuffer(GL_ARRAY_BUFFER,
                         m_colorBufferObject->getBufferObjectName());
            glBufferData(GL_ARRAY_BUFFER,
                         colorSizeBytes,
                         colorDataPointer,
                         usageHint);
        }
            break;
        case GraphicsPrimitive::ColorDataType::UNSIGNED_BYTE_RGBA:
        {
            EventGraphicsOpenGLCreateBufferObject createEvent;
            EventManager::get()->sendEvent(createEvent.getPointer());
            m_colorBufferObject.reset(createEvent.getOpenGLBufferObject());
            CaretAssert(m_colorBufferObject->getBufferObjectName());
            
            m_componentsPerColor = 4;
            m_colorDataType = GL_UNSIGNED_BYTE;
            
            const GLuint colorSizeBytes = primitive->m_unsignedByteRGBA.size() * sizeof(uint8_t);
            CaretAssert(colorSizeBytes > 0);
            const GLvoid* colorDataPointer = (const GLvoid*)&primitive->m_unsignedByteRGBA[0];
            
            glBindBuffer(GL_ARRAY_BUFFER,
                         m_colorBufferObject->getBufferObjectName());
            glBufferData(GL_ARRAY_BUFFER,
                         colorSizeBytes,
                         colorDataPointer,
                         usageHint);
        }
            break;
    }
    
    m_reloadColorsFlag = false;
}

/**
 * Load the texture coordinate buffer.
 * @param primitive
 *     The graphics primitive that will be drawn.
 */
void
GraphicsEngineDataOpenGL::loadTextureCoordinateBuffer(GraphicsPrimitive* primitive)
{
    CaretAssert(primitive);
    
    GLenum usageHint = getOpenGLBufferUsageHint(primitive->getUsageTypeTextureCoordinates());
    
    switch (primitive->m_textureDataType) {
        case GraphicsPrimitive::TextureDataType::FLOAT_STR:
        {
            EventGraphicsOpenGLCreateBufferObject createEvent;
            EventManager::get()->sendEvent(createEvent.getPointer());
            m_textureCoordinatesBufferObject.reset(createEvent.getOpenGLBufferObject());
            CaretAssert(m_textureCoordinatesBufferObject->getBufferObjectName());
            
            m_textureCoordinatesDataType = GL_FLOAT;
            const GLuint textureSizeBytes = primitive->m_floatTextureSTR.size() * sizeof(float);
            CaretAssert(textureSizeBytes > 0);
            const GLvoid* textureDataPointer = (const GLvoid*)&primitive->m_floatTextureSTR[0];
            
            glBindBuffer(GL_ARRAY_BUFFER,
                         m_textureCoordinatesBufferObject->getBufferObjectName());
            glBufferData(GL_ARRAY_BUFFER,
                         textureSizeBytes,
                         textureDataPointer,
                         usageHint);
            
        }
            break;
        case GraphicsPrimitive::TextureDataType::NONE:
            break;
    }    
}


/**
 * Load the buffers with data from the grpahics primitive.
 *
 * @param primitive
 *     The graphics primitive that will be drawn.
 */
void
GraphicsEngineDataOpenGL::loadAllBuffers(GraphicsPrimitive* primitive)
{
    loadCoordinateBuffer(primitive);
    loadNormalVectorBuffer(primitive);
    loadColorBuffer(primitive);
    loadTextureCoordinateBuffer(primitive);    
}

/**
 * If needed, load the texture buffer.
 * Note that the texture buffer may be invalidated when
 * an image is captures as the image capture wipes out all
 * textures.
 *
 * @param primitive
 *     The graphics primitive that will be drawn.
*/
void
GraphicsEngineDataOpenGL::loadTextureImageDataBuffer(GraphicsPrimitive* primitive)
{
    if (m_textureImageDataName != NULL) {
        return;
    }
    
    switch (primitive->m_textureDataType) {
        case GraphicsPrimitive::TextureDataType::FLOAT_STR:
        {
            const int32_t imageWidth  = primitive->m_textureImageWidth;
            const int32_t imageHeight = primitive->m_textureImageHeight;
            const int32_t imageNumberOfBytes = imageWidth * imageHeight * 4;
            if (imageNumberOfBytes <= 0) {
                CaretLogWarning("Invalid texture (empty data) for drawing primitive with texture");
                return;
            }
            if (imageNumberOfBytes != static_cast<int32_t>(primitive->m_textureImageBytesRGBA.size())) {
                CaretLogWarning("Image bytes size incorrect for image width/height");
                return;
            }
            
            glPushClientAttrib(GL_CLIENT_PIXEL_STORE_BIT);
            
            const GLubyte* imageBytesRGBA = &primitive->m_textureImageBytesRGBA[0];
            EventGraphicsOpenGLCreateTextureName createEvent;
            EventManager::get()->sendEvent(createEvent.getPointer());
            m_textureImageDataName = createEvent.getOpenGLTextureName();
            CaretAssert(m_textureImageDataName);
            
            const GLuint openGLTextureName = m_textureImageDataName->getTextureName();
            
            glBindTexture(GL_TEXTURE_2D, openGLTextureName);
            
            bool useMipMapFlag = true;
            switch (primitive->getTextureFilteringType()) {
                case GraphicsPrimitive::TextureFilteringType::LINEAR:
                    break;
                case GraphicsPrimitive::TextureFilteringType::NEAREST:
                    useMipMapFlag = false;
                    break;
            }
            
            if (useMipMapFlag) {
                switch (primitive->getTextureWrappingType()) {
                    case GraphicsPrimitive::TextureWrappingType::CLAMP:
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
                        break;
                    case GraphicsPrimitive::TextureWrappingType::REPEAT:
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
                        break;
                }
                
                switch (primitive->getTextureFilteringType()) {
                    case GraphicsPrimitive::TextureFilteringType::LINEAR:
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR_MIPMAP_LINEAR);
                        break;
                    case GraphicsPrimitive::TextureFilteringType::NEAREST:
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
                        useMipMapFlag = false;
                        break;
                }
                
                /*
                 * This code seems to work if OpenGL 3.0 or later and
                 * replaces gluBuild2DMipmaps()
                 *
                 * glTexImage2D(GL_TEXTURE_2D,     // MUST BE GL_TEXTURE_2D
                 *            0,                 // level of detail 0=base, n is nth mipmap reduction
                 *            GL_RGBA,           // number of components
                 *            imageWidth,        // width of image
                 *            imageHeight,       // height of image
                 *            0,                 // border
                 *            GL_RGBA,           // format of the pixel data
                 *            GL_UNSIGNED_BYTE,  // data type of pixel data
                 *            imageBytesRGBA);   // pointer to image data
                 * glGenerateMipmap(GL_TEXTURE_2D);
                 */
                
                
                const int errorCode = gluBuild2DMipmaps(GL_TEXTURE_2D,     // MUST BE GL_TEXTURE_2D
                                                        GL_RGBA,           // number of components
                                                        imageWidth,        // width of image
                                                        imageHeight,       // height of image
                                                        GL_RGBA,           // format of the pixel data
                                                        GL_UNSIGNED_BYTE,  // data type of pixel data
                                                        imageBytesRGBA);    // pointer to image data
                if (errorCode != 0) {
                    useMipMapFlag = false;
                    
                    const GLubyte* errorChars = gluErrorString(errorCode);
                    if (errorChars != NULL) {
                        const QString errorText = ("ERROR building mipmaps for annotation image: "
                                                   + AString((char*)errorChars));
                        CaretLogSevere(errorText);
                    }
                }
            }
            
            if ( ! useMipMapFlag) {
                switch (primitive->getTextureWrappingType()) {
                    case GraphicsPrimitive::TextureWrappingType::CLAMP:
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP);
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP);
                        break;
                    case GraphicsPrimitive::TextureWrappingType::REPEAT:
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
                        break;
                }
                
                switch (primitive->getTextureFilteringType()) {
                    case GraphicsPrimitive::TextureFilteringType::LINEAR:
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
                        break;
                    case GraphicsPrimitive::TextureFilteringType::NEAREST:
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
                        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
                        break;
                }
                
                glTexImage2D(GL_TEXTURE_2D,     // MUST BE GL_TEXTURE_2D
                             0,                 // level of detail 0=base, n is nth mipmap reduction
                             GL_RGBA,           // number of components
                             imageWidth,        // width of image
                             imageHeight,       // height of image
                             0,                 // border
                             GL_RGBA,           // format of the pixel data
                             GL_UNSIGNED_BYTE,  // data type of pixel data
                             imageBytesRGBA);   // pointer to image data
                
                const auto errorGL = GraphicsUtilitiesOpenGL::getOpenGLError();
                if (errorGL) {
                    CaretLogSevere("OpenGL error after glTexImage2D(), width="
                                   + AString::number(imageWidth)
                                   + ", height="
                                   + AString::number(imageHeight)
                                   + ": "
                                   + errorGL->getVerboseDescription());
                }
            }

            glBindTexture(GL_TEXTURE_2D, 0);
            glPopClientAttrib();
        }
            break;
        case GraphicsPrimitive::TextureDataType::NONE:
            break;
    }
}

/**
 * Draw the given graphics primitive.
 *
 * @param primitive
 *     Primitive that is drawn.
 */
void
GraphicsEngineDataOpenGL::draw(GraphicsPrimitive* primitive)
{
    CaretAssert(primitive);
    
    if ( ! primitive->isValid()) {
        CaretLogFine("Attempting to draw invalid Graphics Primitive");
        return;
    }
    
    GraphicsPrimitive* primitiveToDraw = primitive;
    
    /* Conversion to window space creates a new primitive that must be deleted */
    std::unique_ptr<GraphicsPrimitive> lineConversionPrimitive;
    
    bool modelSpaceLineFlag = false;
    bool windowSpaceLineFlag = false;
    bool millimeterPointsFlag = false;
    bool spheresFlag = false;
    switch (primitive->m_primitiveType) {
        case GraphicsPrimitive::PrimitiveType::OPENGL_LINE_LOOP:
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_LINE_STRIP:
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_LINES:
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_POINTS:
            switch (primitive->m_pointSizeType) {
                case GraphicsPrimitive::PointSizeType::MILLIMETERS:
                    millimeterPointsFlag = true;
                    break;
                case GraphicsPrimitive::PointSizeType::PERCENTAGE_VIEWPORT_HEIGHT:
                    break;
                case GraphicsPrimitive::PointSizeType::PIXELS:
                    break;
            }
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_TRIANGLE_FAN:
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_TRIANGLE_STRIP:
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_TRIANGLES:
            break;
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_MITER_JOIN:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_MITER_JOIN:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINES:
            modelSpaceLineFlag = true;
            break;
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_LOOP_MITER_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_STRIP_MITER_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINES:
            windowSpaceLineFlag = true;
            break;
        case GraphicsPrimitive::PrimitiveType::SPHERES:
            spheresFlag = true;
            break;
    }
    
    if (modelSpaceLineFlag
        || windowSpaceLineFlag) {
        if (primitive->getNumberOfVertices() == 1) {
            spheresFlag = true;
            
            uint8_t rgba[4];
            switch (primitive->getColorDataType()) {
                case GraphicsPrimitive::ColorDataType::FLOAT_RGBA:
                {
                    float rgbaFloat[4];
                    primitive->getVertexFloatRGBA(0, rgbaFloat);
                    for (int32_t i = 0; i < 4; i++) {
                        rgba[i] = static_cast<uint8_t>(rgbaFloat[i] / 255.0);
                    }
                }
                    break;
                case GraphicsPrimitive::ColorDataType::NONE:
                    break;
                case GraphicsPrimitive::ColorDataType::UNSIGNED_BYTE_RGBA:
                    primitive->getVertexByteRGBA(0, rgba);
                    break;
            }
            float sphereXYZ[3];
            primitive->getVertexFloatXYZ(0, sphereXYZ);
            GraphicsPrimitive::LineWidthType lineWidthType;
            float lineWidth;
            primitive->getLineWidth(lineWidthType, lineWidth);
            std::unique_ptr<GraphicsPrimitiveV3f> pointPrimitive(GraphicsPrimitive::newPrimitiveV3f(GraphicsPrimitive::PrimitiveType::OPENGL_POINTS,
                                                                                                     rgba));
            pointPrimitive->addVertex(sphereXYZ);
            
            GraphicsPrimitive::PointSizeType pointSizeType(GraphicsPrimitive::PointSizeType::PIXELS);
            switch (lineWidthType) {
                case GraphicsPrimitive::LineWidthType::PERCENTAGE_VIEWPORT_HEIGHT:
                    pointSizeType = GraphicsPrimitive::PointSizeType::PERCENTAGE_VIEWPORT_HEIGHT;
                    break;
                case GraphicsPrimitive::LineWidthType::PIXELS:
                    pointSizeType = GraphicsPrimitive::PointSizeType::PIXELS;
                    break;
            }
            pointPrimitive->setPointDiameter(pointSizeType,
                                             lineWidth);

            drawPrivate(PrivateDrawMode::DRAW_NORMAL,
                        pointPrimitive.get(),
                        NULL);
            return;
        }
    }
    
    if (millimeterPointsFlag) {
        /* 
         * Special case for points in millimeters 
         */
        drawPointsPrimitiveMillimeters(primitive);
    }
    else if (spheresFlag) {
        drawSpheresPrimitive(primitive);
    }
    else if (modelSpaceLineFlag
             || windowSpaceLineFlag) {
        AString errorMessage;
        primitiveToDraw = GraphicsOpenGLPolylineTriangles::convertWorkbenchLinePrimitiveTypeToOpenGL(primitive,
                                                                                                 errorMessage);
        if (primitiveToDraw == NULL) {
            const AString msg("For developer: "
                              + errorMessage);
#ifdef NDEBUG
            CaretLogFine(msg);
#else
            CaretLogSevere(msg);
#endif
            return;
        }
        SpaceMode spaceMode = SpaceMode::WINDOW;
        if (modelSpaceLineFlag) {
            spaceMode = SpaceMode::MODEL;
        }
        else if (windowSpaceLineFlag) {
            spaceMode = SpaceMode::WINDOW;
        }
        else {
            CaretAssert(0);
        }
        
        lineConversionPrimitive.reset(primitiveToDraw);
        drawModelOrWindowSpace(spaceMode,
                               PrivateDrawMode::DRAW_NORMAL,
                               primitiveToDraw,
                               NULL);
    }
    else {
        drawPrivate(PrivateDrawMode::DRAW_NORMAL,
                    primitiveToDraw,
                    NULL);
    }
}

/**
 * Draw the graphics primitive in model space.
 *
 * @param spaceMode
 *     Space mode: model or window
 * @param drawMode
 *     Mode for drawing.
 * @param primitive
 *     Primitive that is drawn.
 * @param primitiveSelectionHelper
 *     Selection helper when draw mode is selection.
 */
void
GraphicsEngineDataOpenGL::drawModelOrWindowSpace(const SpaceMode spaceMode,
                                                 const PrivateDrawMode drawMode,
                                                 GraphicsPrimitive* primitive,
                                                 GraphicsPrimitiveSelectionHelper* primitiveSelectionHelper)
{
    bool windowSpaceFlag = false;
    switch (spaceMode) {
        case MODEL:
            break;
        case WINDOW:
            windowSpaceFlag = true;
            break;
    }
    
    int32_t polygonMode[2];
    int32_t viewport[4];
    saveOpenGLStateForWindowSpaceDrawing(polygonMode,
                                         viewport);
    
    if (windowSpaceFlag) {
        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        glOrtho(viewport[0], viewport[0] + viewport[2],
                viewport[1], viewport[1] + viewport[3],
                0, 1);
        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
    }
    
    glPolygonMode(GL_FRONT, GL_FILL);
    
    drawPrivate(drawMode,
                primitive,
                primitiveSelectionHelper);
    
    restoreOpenGLStateForWindowSpaceDrawing(polygonMode,
                                            viewport);
}

/**
 * Draw a point primitive type in millimeters (same as model space).
 *
 * @param primitive
 *     The points primitive type.
 */
void
GraphicsEngineDataOpenGL::drawPointsPrimitiveMillimeters(const GraphicsPrimitive* primitive)
{
    CaretAssert(primitive);
    CaretAssert(primitive->getPrimitiveType() == GraphicsPrimitive::PrimitiveType::OPENGL_POINTS);
    CaretAssert(primitive->m_pointSizeType == GraphicsPrimitive::PointSizeType::MILLIMETERS);
    
    GraphicsPrimitive::PointSizeType sizeType;
    float sizeValue = 0.0f;
    primitive->getPointDiameter(sizeType,
                                sizeValue);
    
    const int32_t numberOfVertices = primitive->getNumberOfVertices();
    switch (primitive->m_pointSizeType) {
        case GraphicsPrimitive::PointSizeType::MILLIMETERS:
        {
            switch (primitive->m_vertexColorType) {
                case GraphicsPrimitive::VertexColorType::NONE:
                    CaretAssert(0);
                    CaretLogSevere("NONE vertex color type, should never occur when drawing points in millimeters");
                    return;
                    break;
                case GraphicsPrimitive::VertexColorType::PER_VERTEX_RGBA:
                {
                    for (int32_t i = 0; i < numberOfVertices; i++) {
                        const int32_t i3 = i * 3;
                        CaretAssertVectorIndex(primitive->m_xyz, i3 + 2);
                        const float* xyz = &primitive->m_xyz[i3];
                        const int32_t i4 = i * 4;
                        
                        uint8_t* rgba = NULL;
                        switch (primitive->m_colorDataType) {
                            case GraphicsPrimitive::ColorDataType::FLOAT_RGBA:
                                CaretAssert(0);
                                break;
                            case GraphicsPrimitive::ColorDataType::UNSIGNED_BYTE_RGBA:
                                CaretAssertVectorIndex(primitive->m_unsignedByteRGBA, i4 + 3);
                                rgba = const_cast<uint8_t*>(&primitive->m_unsignedByteRGBA[i4]);
                                break;
                            case GraphicsPrimitive::ColorDataType::NONE:
                                CaretAssert(0);
                                break;
                        }
                        
                        GraphicsShape::drawSquare(xyz, rgba, sizeValue);
                    }
                }
                    break;
                case GraphicsPrimitive::VertexColorType::SOLID_RGBA:
                {
                    const int32_t numberOfVertices = primitive->getNumberOfVertices();
                    CaretAssertVectorIndex(primitive->m_xyz, (numberOfVertices - 1) * 3 + 1);
                    GraphicsShape::drawSquares(&primitive->m_xyz[0],
                                               numberOfVertices,
                                               &primitive->m_unsignedByteRGBA[0],
                                               sizeValue);
                }
                    break;
            }
        }
            break;
        case GraphicsPrimitive::PointSizeType::PERCENTAGE_VIEWPORT_HEIGHT:
        case GraphicsPrimitive::PointSizeType::PIXELS:
            CaretAssert(0);
            break;
    }
}

/**
 * Draw a sphere primitive type.
 * 
 * @param primitive
 *     The spheres primitive type.
 */
void
GraphicsEngineDataOpenGL::drawSpheresPrimitive(const GraphicsPrimitive* primitive)
{
    CaretAssert(primitive);
    CaretAssert(primitive->getPrimitiveType() == GraphicsPrimitive::PrimitiveType::SPHERES);

    GraphicsPrimitive::SphereSizeType sizeType;
    float sizeValue = 0.0f;
    primitive->getSphereDiameter(sizeType,
                                 sizeValue);
    
    switch (primitive->m_vertexColorType) {
        case GraphicsPrimitive::VertexColorType::NONE:
            CaretAssert(0);
            CaretLogSevere("NONE vertex color type, should never occur when drawing spheres");
            return;
            break;
        case GraphicsPrimitive::VertexColorType::PER_VERTEX_RGBA:
        {
            const int32_t numberOfVertices = primitive->getNumberOfVertices();
            for (int32_t i = 0; i < numberOfVertices; i++) {
                const int32_t i3 = i * 3;
                CaretAssertVectorIndex(primitive->m_xyz, i3 + 2);
                const int32_t i4 = i * 4;
                
                switch (primitive->m_colorDataType) {
                    case GraphicsPrimitive::ColorDataType::FLOAT_RGBA:
                        CaretAssert(0);
                        break;
                    case GraphicsPrimitive::ColorDataType::UNSIGNED_BYTE_RGBA:
                        CaretAssertVectorIndex(primitive->m_unsignedByteRGBA, i4 + 3);
                        GraphicsShape::drawSphereByteColor(&primitive->m_xyz[i3],
                                                           &primitive->m_unsignedByteRGBA[i4],
                                                           sizeValue);
                        break;
                    case GraphicsPrimitive::ColorDataType::NONE:
                        CaretAssert(0);
                        break;
                }
            }
        }
            break;
        case GraphicsPrimitive::VertexColorType::SOLID_RGBA:
        {
            const int32_t numberOfVertices = primitive->getNumberOfVertices();
            CaretAssertVectorIndex(primitive->m_xyz, (numberOfVertices - 1) * 3 + 1);
            GraphicsShape::drawSpheresByteColor(&primitive->m_xyz[0],
                                                numberOfVertices,
                                                &primitive->m_unsignedByteRGBA[0],
                                                sizeValue);
        }
            break;
    }
}

/**
 * Save the OpenGL State for drawing in window space.
 *
 * @param polygonMode
 *    Output with the polygon mode.
 * @param viewport
 *    Output with viewport.
 */
void
GraphicsEngineDataOpenGL::saveOpenGLStateForWindowSpaceDrawing(int32_t polygonMode[2],
                                          int32_t viewport[4])
{
    glPushAttrib(GL_COLOR_BUFFER_BIT
                 | GL_ENABLE_BIT
                 | GL_HINT_BIT
                 | GL_LIGHTING_BIT
                 | GL_POLYGON_BIT);
    
    /*
     * For anti-aliasing requires sorting polygons
     * and turning of depth testing
     */
    const bool drawWithAntiAliasingFlag = false;
    if (drawWithAntiAliasingFlag) {
        glBlendFunc (GL_SRC_ALPHA_SATURATE, GL_ONE);
        glEnable(GL_BLEND);
        glHint(GL_POLYGON_SMOOTH_HINT,
               GL_NICEST);
        glEnable(GL_POLYGON_SMOOTH);
    }
    
    glDisable(GL_LIGHTING);
    glDisable(GL_COLOR_MATERIAL);
    glGetIntegerv(GL_POLYGON_MODE,
                  polygonMode);
    
    glGetIntegerv(GL_VIEWPORT,
                  viewport);
    glMatrixMode(GL_PROJECTION);
    glPushMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPushMatrix();
}

/**
 * Restore the OpenGL State from drawing in window space.
 
 * @param polygonMode
 *    Polygon mode.
 * @param viewport
 *    Viewport.
 */
void
GraphicsEngineDataOpenGL::restoreOpenGLStateForWindowSpaceDrawing(int32_t polygonMode[2],
                                                                  int32_t viewport[4])
{
    glPolygonMode(GL_FRONT, polygonMode[0]);
    glPolygonMode(GL_BACK, polygonMode[1]);
    glMatrixMode(GL_PROJECTION);
    glPopMatrix();
    glMatrixMode(GL_MODELVIEW);
    glPopMatrix();
    glViewport(viewport[0], viewport[1], viewport[2], viewport[3]);
    
    glPopAttrib();
}

/**
 * Draw to determine the index of primitive at the given (X, Y) location
 *
 * @param primitive
 *     Primitive that is drawn.
 * @param pixelX
 *     X location for selection.
 * @param pixelY
 *     Y location for selection.
 * @param selectedPrimitiveIndexOut
 *     Output with the index of the primitive (invalid no primitve at location
 * @param selectedPrimitiveDepthOut
 *     Output with depth at the location.
 */
void
GraphicsEngineDataOpenGL::drawWithSelection(GraphicsPrimitive* primitive,
                                            const int32_t pixelX,
                                            const int32_t pixelY,
                                            int32_t& selectedPrimitiveIndexOut,
                                            float&   selectedPrimitiveDepthOut)
{
    CaretAssert(primitive);
    
    bool modelSpaceLineFlag = false;
    bool windowSpaceLineFlag = false;
    switch (primitive->m_primitiveType) {
        case GraphicsPrimitive::PrimitiveType::OPENGL_LINE_LOOP:
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_LINE_STRIP:
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_LINES:
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_POINTS:
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_TRIANGLE_FAN:
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_TRIANGLE_STRIP:
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_TRIANGLES:
            break;
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_MITER_JOIN:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_MITER_JOIN:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINES:
            modelSpaceLineFlag = true;
            break;
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_LOOP_MITER_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_STRIP_MITER_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINES:
            windowSpaceLineFlag = true;
            break;
        case GraphicsPrimitive::PrimitiveType::SPHERES:
            CaretAssertMessage(0, "Not yet implemented");
            break;
    }
    
    selectedPrimitiveIndexOut = -1;
    selectedPrimitiveDepthOut = 0.0;

    GLint viewport[4];
    glGetIntegerv(GL_VIEWPORT,
                  viewport);
    
    /*
     * Must clear color and depth buffers
     */
    glClear(GL_COLOR_BUFFER_BIT
            | GL_DEPTH_BUFFER_BIT);
    
    if ((pixelX    >= viewport[0])
        && (pixelX <  (viewport[0] + viewport[2]))
        && (pixelY >= viewport[1])
        && (pixelY <  (viewport[1] + viewport[3]))) {
        /*
         * Saves glPixelStore parameters
         */
        glPushClientAttrib(GL_CLIENT_PIXEL_STORE_BIT);
        glPushAttrib(GL_LIGHTING_BIT);
        
        glShadeModel(GL_FLAT);
        glDisable(GL_LIGHTING);
        
        std::vector<int32_t> triangleVertexIndicesToLineVertexIndices;
        std::unique_ptr<GraphicsPrimitiveSelectionHelper> selectionHelper;
        if (modelSpaceLineFlag
            || windowSpaceLineFlag) {
            AString errorMessage;
            
            std::unique_ptr<GraphicsPrimitive> primitiveToDraw(GraphicsOpenGLPolylineTriangles::convertWorkbenchLinePrimitiveTypeToOpenGL(primitive,
                                                                                                                                          &triangleVertexIndicesToLineVertexIndices,
                                                                                                                                          errorMessage));
            if (primitiveToDraw == NULL) {
                const AString msg("For developer: "
                                  + errorMessage);
#ifdef NDEBUG
                CaretLogFine(msg);
#else
                CaretLogSevere(msg);
#endif
                return;
            }
            SpaceMode spaceMode = SpaceMode::WINDOW;
            if (modelSpaceLineFlag) {
                spaceMode = SpaceMode::MODEL;
            }
            else if (windowSpaceLineFlag) {
                spaceMode = SpaceMode::WINDOW;
            }
            else {
                CaretAssert(0);
            }
            
            selectionHelper.reset(new GraphicsPrimitiveSelectionHelper(primitiveToDraw.get()));
            selectionHelper->setupSelectionBeforeDrawing();
            drawModelOrWindowSpace(spaceMode,
                                   PrivateDrawMode::DRAW_SELECTION,
                                   primitiveToDraw.get(),
                                   selectionHelper.get());
        }
        else {
            selectionHelper.reset(new GraphicsPrimitiveSelectionHelper(primitive));
            selectionHelper->setupSelectionBeforeDrawing();
            drawPrivate(PrivateDrawMode::DRAW_SELECTION,
                        primitive,
                        selectionHelper.get());
        }
        
        /*
         * QOpenGLWidget Note: The QOpenGLWidget always renders in a
         * frame buffer object (see its documentation).  This is
         * probably why calls to glReadBuffer() always cause an
         * OpenGL error.
         */
#ifdef WORKBENCH_USE_QT5_QOPENGL_WIDGET
        /* do not call glReadBuffer() */
#else
        glReadBuffer(GL_BACK);
#endif
        uint8_t pixelRGBA[4];

        glPixelStorei(GL_PACK_SKIP_ROWS, 0);
        glPixelStorei(GL_PACK_SKIP_PIXELS, 0);
        glPixelStorei(GL_PACK_ALIGNMENT, 1); /* bytes */
        glReadPixels(pixelX,
                     pixelY,
                     1,
                     1,
                     GL_RGBA,
                     GL_UNSIGNED_BYTE,
                     pixelRGBA);
        
        /*
         * Get depth from depth buffer
         */
        glPixelStorei(GL_PACK_ALIGNMENT, 4); /* float */
        glReadPixels(pixelX,
                     pixelY,
                     1,
                     1,
                     GL_DEPTH_COMPONENT,
                     GL_FLOAT,
                     &selectedPrimitiveDepthOut);
        
        glPopAttrib();
        glPopClientAttrib();
        
        CaretAssert(selectionHelper.get());
        selectedPrimitiveIndexOut = selectionHelper->getPrimitiveIndexFromEncodedRGBA(pixelRGBA);
        
        if (modelSpaceLineFlag
            || windowSpaceLineFlag) {
            if ((selectedPrimitiveIndexOut >= 0)
                && selectedPrimitiveIndexOut < static_cast<int32_t>(triangleVertexIndicesToLineVertexIndices.size())) {
                CaretAssertVectorIndex(triangleVertexIndicesToLineVertexIndices, selectedPrimitiveIndexOut);
                selectedPrimitiveIndexOut = triangleVertexIndicesToLineVertexIndices[selectedPrimitiveIndexOut];
                
                if (selectedPrimitiveIndexOut >= primitive->getNumberOfVertices()) {
                    selectedPrimitiveIndexOut = -1;
                }
            }
            else {
                selectedPrimitiveIndexOut = -1;
            }
        }
        else {
            if (selectedPrimitiveIndexOut >= 0) {
                if (selectedPrimitiveIndexOut >= primitive->getNumberOfVertices()) {
                    selectedPrimitiveIndexOut = -1;
                }
            }
        }
    }
}

/**
 * Draw the graphics primitive.
 *
 * @param drawMode
 *     Mode for drawing.
 * @param primitive
 *     Primitive that is drawn.
 * @param primitiveSelectionHelper
 *     Selection helper when draw mode is selection.
 */
void
GraphicsEngineDataOpenGL::drawPrivate(const PrivateDrawMode drawMode,
                                      GraphicsPrimitive* primitive,
                                      GraphicsPrimitiveSelectionHelper* primitiveSelectionHelper)
{
    CaretAssert(primitive);
    if ( ! primitive->isValid()) {
        return;
    }
    
    switch (drawMode) {
        case PrivateDrawMode::DRAW_NORMAL:
            break;
        case PrivateDrawMode::DRAW_SELECTION:
            break;
    }
    
    GraphicsEngineDataOpenGL* openglData = primitive->getGraphicsEngineDataForOpenGL();
    if (openglData == NULL) {
        openglData = new GraphicsEngineDataOpenGL();
        primitive->setGraphicsEngineDataForOpenGL(openglData);
        
        openglData->loadAllBuffers(primitive);
    }
    else {
        /*
         * Coordinates may get updated.
         */
        if (openglData->m_reloadCoordinatesFlag) {
            openglData->loadCoordinateBuffer(primitive);
        }
        
        /*
         * Colors may get updated
         */
        if (openglData->m_reloadColorsFlag) {
            openglData->loadColorBuffer(primitive);
        }
    }
    
    openglData->loadTextureImageDataBuffer(primitive);
    
    /*
     * Primitive may want to delete its instance data to save memory
     */
    primitive->setOpenGLBuffersHaveBeenLoadedByGraphicsEngine();
    
    const GLuint coordBufferID = openglData->m_coordinateBufferObject->getBufferObjectName();
    CaretAssert(coordBufferID);
    if ( ! glIsBuffer(coordBufferID)) {
        CaretAssertMessage(0, "Coordinate buffer is INVALID");
        CaretLogSevere("Coordinate buffer is INVALID");
        return;
    }
    
    switch (drawMode) {
        case PrivateDrawMode::DRAW_NORMAL:
            break;
        case PrivateDrawMode::DRAW_SELECTION:
            break;
    }
    
    /*
     * Setup vertices for drawing
     */
    glEnableClientState(GL_VERTEX_ARRAY);
    CaretAssert(glIsBuffer(coordBufferID));
    glBindBuffer(GL_ARRAY_BUFFER,
                 coordBufferID);
    glVertexPointer(openglData->m_coordinatesPerVertex,
                    openglData->m_coordinateDataType,
                    0,
                    (GLvoid*)0);
    
    /*
     * Setup normals for drawing
     */
    if (openglData->m_normalVectorBufferObject != NULL) {
        glEnableClientState(GL_NORMAL_ARRAY);
        CaretAssert(glIsBuffer(openglData->m_normalVectorBufferObject->getBufferObjectName()));
        glBindBuffer(GL_ARRAY_BUFFER, openglData->m_normalVectorBufferObject->getBufferObjectName());
        glNormalPointer(openglData->m_normalVectorDataType, 0, (GLvoid*)0);
    }
    else {
        glDisableClientState(GL_NORMAL_ARRAY);
    }
    
    bool hasColorFlag = false;
    switch (primitive->getColorDataType()) {
        case GraphicsPrimitive::ColorDataType::NONE:
            break;
        case GraphicsPrimitive::ColorDataType::FLOAT_RGBA:
            hasColorFlag = true;
            break;
        case GraphicsPrimitive::ColorDataType::UNSIGNED_BYTE_RGBA:
            hasColorFlag = true;
            break;
    }
    /*
     * Setup colors for drawing
     */
    GLuint localColorBufferName = 0;
    if (hasColorFlag) {
        switch (drawMode) {
            case PrivateDrawMode::DRAW_NORMAL:
            {
                glEnableClientState(GL_COLOR_ARRAY);
                CaretAssert(glIsBuffer(openglData->m_colorBufferObject->getBufferObjectName()));
                glBindBuffer(GL_ARRAY_BUFFER, openglData->m_colorBufferObject->getBufferObjectName());
                glColorPointer(openglData->m_componentsPerColor, openglData->m_colorDataType, 0, (GLvoid*)0);
            }
                break;
            case PrivateDrawMode::DRAW_SELECTION:
            {
                const int32_t numberOfVertices = primitive->m_xyz.size() / 3;
                if (numberOfVertices > 0) {
                    glGenBuffers(1, &localColorBufferName);
                    CaretAssert(localColorBufferName > 0);
                    
                    const GLint  componentsPerColor = 4;
                    const GLenum colorDataType = GL_UNSIGNED_BYTE;
                    
                    const std::vector<uint8_t> selectionRGBA = primitiveSelectionHelper->getSelectionEncodedRGBA();
                    if (static_cast<int32_t>(selectionRGBA.size() / 4) == numberOfVertices) {
                        const GLuint colorSizeBytes = selectionRGBA.size() * sizeof(GLubyte);
                        const GLvoid* colorDataPointer = (const GLvoid*)&selectionRGBA[0];
                        glBindBuffer(GL_ARRAY_BUFFER,
                                     localColorBufferName);
                        glBufferData(GL_ARRAY_BUFFER,
                                     colorSizeBytes,
                                     colorDataPointer,
                                     GL_STREAM_DRAW);
                        glEnableClientState(GL_COLOR_ARRAY);
                        glColorPointer(componentsPerColor, colorDataType, 0, (GLvoid*)0);
                    }
                }
            }
                break;
        }
    }
    
    bool hasTextureFlag = false;
    switch (primitive->getTextureDataType()) {
        case GraphicsPrimitive::TextureDataType::NONE:
            break;
        case GraphicsPrimitive::TextureDataType::FLOAT_STR:
            hasTextureFlag = true;
            break;
    }
    
    /*
     * Setup textures for drawing
     */
    if (hasTextureFlag
        && (drawMode == PrivateDrawMode::DRAW_NORMAL)
        && (openglData->m_textureImageDataName->getTextureName() > 0)) {
            glEnable(GL_TEXTURE_2D);
            glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
            glBindTexture(GL_TEXTURE_2D, openglData->m_textureImageDataName->getTextureName());
            
            glEnableClientState(GL_TEXTURE_COORD_ARRAY);
            CaretAssert(glIsBuffer(openglData->m_textureCoordinatesBufferObject->getBufferObjectName()));
            glBindBuffer(GL_ARRAY_BUFFER, openglData->m_textureCoordinatesBufferObject->getBufferObjectName());
        glTexCoordPointer(3, openglData->m_textureCoordinatesDataType, 0, (GLvoid*)0);
    }
     
    /*
     * Unbind buffers
     */
    glBindBuffer(GL_ARRAY_BUFFER,
                 0);
    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER,
                 0);
    
    /*
     * Draw the primitive
     */
    GLenum openGLPrimitiveType = GL_INVALID_ENUM;
    switch (primitive->m_primitiveType) {
        case GraphicsPrimitive::PrimitiveType::OPENGL_LINE_LOOP:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_LOOP_MITER_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_LOOP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_LOOP_MITER_JOIN:
            openGLPrimitiveType = GL_LINE_LOOP;
            glLineWidth(getLineWidthForDrawingInPixels(primitive));
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_LINE_STRIP:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINE_STRIP_MITER_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_STRIP_BEVEL_JOIN:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINE_STRIP_MITER_JOIN:
            openGLPrimitiveType = GL_LINE_STRIP;
            glLineWidth(getLineWidthForDrawingInPixels(primitive));
            break;
        case GraphicsPrimitive::PrimitiveType::MODEL_SPACE_POLYGONAL_LINES:
        case GraphicsPrimitive::PrimitiveType::OPENGL_LINES:
        case GraphicsPrimitive::PrimitiveType::POLYGONAL_LINES:
            openGLPrimitiveType = GL_LINES;
            glLineWidth(getLineWidthForDrawingInPixels(primitive));
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_POINTS:
            openGLPrimitiveType = GL_POINTS;
            glPointSize(getPointDiameterForDrawingInPixels(primitive));
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_TRIANGLE_FAN:
            openGLPrimitiveType = GL_TRIANGLE_FAN;
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_TRIANGLE_STRIP:
            openGLPrimitiveType = GL_TRIANGLE_STRIP;
            break;
        case GraphicsPrimitive::PrimitiveType::OPENGL_TRIANGLES:
            openGLPrimitiveType = GL_TRIANGLES;
            break;
        case GraphicsPrimitive::PrimitiveType::SPHERES:
            CaretAssertMessage(0, "SPHERES are not drawn with this method");
            break;
    }
    
    CaretAssert(openGLPrimitiveType != GL_INVALID_ENUM);
    glDrawArrays(openGLPrimitiveType, 0, openglData->m_arrayIndicesCount);
    
    /*
     * Disable drawing
     */
    glDisableClientState(GL_VERTEX_ARRAY);
    glDisableClientState(GL_NORMAL_ARRAY);
    glDisableClientState(GL_COLOR_ARRAY);
    glDisableClientState(GL_TEXTURE_COORD_ARRAY);
    
    if (localColorBufferName > 0) {
        glDeleteBuffers(1, &localColorBufferName);
    }
    
    glBindTexture(GL_TEXTURE_2D, 0);
    glDisable(GL_TEXTURE_2D);
}

/**
 * Get the line width in pixels for drawing lines.
 * If the primitive's line width type is a percentage
 * of viewport height, convert it to pixels.
 *
 * @param primtive
 *     The graphics primitive.
 * @return
 *     Width of line in pixels.
 */
float
GraphicsEngineDataOpenGL::getLineWidthForDrawingInPixels(const GraphicsPrimitive* primitive)
{
    CaretAssert(primitive);
    
    float width = primitive->m_lineWidthValue;
    
    switch (primitive->m_lineWidthType) {
        case GraphicsPrimitive::LineWidthType::PERCENTAGE_VIEWPORT_HEIGHT:
        {
            GLint viewport[4];
            glGetIntegerv(GL_VIEWPORT, viewport);
            
            width = (viewport[3] * (width / 100.0f));
        }
            break;
        case GraphicsPrimitive::LineWidthType::PIXELS:
            break;
    }
    
    return width;
}

/**
 * Get the point diameter in pixels for drawing points.
 * If the primitive's point diameter type is a percentage
 * of viewport height, convert it to pixels.
 *
 * @param primtive
 *     The graphics primitive.
 * @return
 *     Diameter of point in pixels.
 */
float
GraphicsEngineDataOpenGL::getPointDiameterForDrawingInPixels(const GraphicsPrimitive* primitive)
{
    CaretAssert(primitive);
    
    float pointSize = primitive->m_pointDiameterValue;
    
    switch (primitive->m_pointSizeType) {
        case GraphicsPrimitive::PointSizeType::MILLIMETERS:
            CaretAssert(0);  /* millimeters not convertible to pixels */
            break;
        case GraphicsPrimitive::PointSizeType::PERCENTAGE_VIEWPORT_HEIGHT:
        {
            GLint viewport[4];
            glGetIntegerv(GL_VIEWPORT, viewport);
            
            pointSize = (viewport[3] * (pointSize / 100.0f));
        }
            break;
        case GraphicsPrimitive::PointSizeType::PIXELS:
            break;
    }
    
    return pointSize;
}

