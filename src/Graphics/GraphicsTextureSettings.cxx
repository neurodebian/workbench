
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

#define __GRAPHICS_TEXTURE_SETTINGS_DECLARE__
#include "GraphicsTextureSettings.h"
#undef __GRAPHICS_TEXTURE_SETTINGS_DECLARE__

#include <cstdint>

#include "CaretAssert.h"
#include "CaretLogger.h"

using namespace caret;


    
/**
 * \class caret::GraphicsTextureSettings 
 * \brief Contains settins used for primitives with textures
 * \ingroup Graphics
 */

/**
 * Static method for allocating a shared pointer pointing to storage for image's RGBA coloring
 * @param imageWidth
 *    Width of image in pixels
 * @param imageHeight
 *    Height of image in pixels
 * @param imageSlices
 *    Number of slices in image (1 if 2D)
 * @param optionalNumberOfBytesOut
 *    If not NULL, contains the number of bytes allocated upon exit
 * @return
 *    Shared point containing memory for loading of rgba color components
 */
std::shared_ptr<uint8_t>
GraphicsTextureSettings::allocateImageRgbaData(const int64_t imageWidth,
                                               const int64_t imageHeight,
                                               const int64_t imageSlices,
                                               int64_t* optionalNumberOfBytesOut)
{
    CaretAssert(imageWidth > 0);
    CaretAssert(imageHeight > 0);
    CaretAssert(imageSlices > 0);
    
    const int64_t bytesPerPixel(4); /* RGBA components*/
    const int64_t numTextureBytes(imageWidth * imageHeight * imageSlices * bytesPerPixel);

    if (optionalNumberOfBytesOut != NULL) {
        *optionalNumberOfBytesOut = numTextureBytes;
    }
    
    /*
     * Compilers before C++17 need deleter for array
     * https://stackoverflow.com/questions/13061979/shared-ptr-to-an-array-should-it-be-used
     * https://kezunlin.me/post/b82753fc/
     * std::shared_ptr<int> sp(new int[10], std::default_delete<int[]>());
     */
    std::shared_ptr<uint8_t> imageData(new uint8_t[numTextureBytes], std::default_delete<uint8_t[]>());
    return imageData;
}

/**
 * Default constructor
 */
GraphicsTextureSettings::GraphicsTextureSettings()
: CaretObject(),
m_imageDataType(ImageDataType::INVALID),
m_imageBytesPointer(NULL)
{
    m_borderColor.fill(0.0f);
}

/**
 * Constructor.
 * @param imageRgbaData
 *    Shared pointer pointing to image data rgba coloring
 *    Use GraphicsTextureSettings::allocateImageRgbaData() to allocate this memory CORRECTLY
 * @param imageWidth
 *    Width of image in pixels
 * @param imageHeight
 *    Height of image in pixels
 * @param imageSlices
 *    Number of slices in image (1 if 2D)
 * @param dimensionType
 *    Texture dimension (2D or 3D)
 * @param pixelFormatType
 *    Pixel color components type
 * @param pixelOrigin
 *    Origin of first pixel in image
 * @param wrappingType
 *    Texture wrapping type
 * @param mipMappingType
 *    Mip mapping type
 * @param compressionType
 *    Type of compression
 * @param magnificationFilter
 *    Type of texture magnification
 * @param minificationFilter
 *    Type of texture minification
 * @param borderColor
 *    Color for texture border (used where there is no texture data)
 */
GraphicsTextureSettings::GraphicsTextureSettings(std::shared_ptr<uint8_t>& imageRgbaData,
                                                 const int64_t         imageWidth,
                                                 const int64_t         imageHeight,
                                                 const int64_t         imageSlices,
                                                 const DimensionType   dimensionType,
                                                 const PixelFormatType pixelFormatType,
                                                 const PixelOrigin     pixelOrigin,
                                                 const WrappingType    wrappingType,
                                                 const MipMappingType  mipMappingType,
                                                 const CompressionType compressionType,
                                                 const GraphicsTextureMagnificationFilterEnum::Enum magnificationFilter,
                                                 const GraphicsTextureMinificationFilterEnum::Enum minificationFilter,
                                                 const std::array<float, 4>& borderColor)
: CaretObject(),
m_imageDataType(ImageDataType::SHARED_PTR),
m_imageRgbaData(imageRgbaData),
m_imageBytesPointer(NULL),
m_imageWidth(imageWidth),
m_imageHeight(imageHeight),
m_imageSlices(imageSlices),
m_dimensionType(dimensionType),
m_pixelFormatType(pixelFormatType),
m_pixelOrigin(pixelOrigin),
m_wrappingType(wrappingType),
m_mipMappingType(mipMappingType),
m_compressionType(compressionType),
m_magnificationFilter(magnificationFilter),
m_minificationFilter(minificationFilter),
m_borderColor(borderColor)
{
    switch (m_mipMappingType) {
        case MipMappingType::DISABLED:
        {
            AString filterName;
            switch (m_minificationFilter) {
                case GraphicsTextureMinificationFilterEnum::LINEAR:
                    break;
                case GraphicsTextureMinificationFilterEnum::LINEAR_MIPMAP_LINEAR:
                    filterName = "LINEAR_MIPMAP_LINEAR";
                    break;
                case GraphicsTextureMinificationFilterEnum::LINEAR_MIPMAP_NEAREST:
                    filterName = "LINEAR_MIPMAP_NEAREST";
                    break;
                case GraphicsTextureMinificationFilterEnum::NEAREST:
                    break;
                case GraphicsTextureMinificationFilterEnum::NEAREST_MIPMAP_LINEAR:
                    filterName = "NEAREST_MIPMAP_LINEAR";
                    break;
                case GraphicsTextureMinificationFilterEnum::NEAREST_MIPMAP_NEAREST:
                    filterName = "NEAREST_MIPMAP_NEAREST";
                    break;
            }
            
            if ( ! filterName.isEmpty()) {
                const AString msg("Mip mapping is DISABLED but the minification filter is set to a "
                                  "mip mapping type "
                                  + filterName
                                  + ".  This may result in nothing or just white displayed.  Either "
                                  "enable mip mapping or use a non-mip mapping minification filter.");
                CaretLogSevere(msg);
            }
        }
            break;
        case MipMappingType::ENABLED:
            break;
    }
}

/**
 * Constructor.
 * @param imageBytesPointer
 *    Pointer to image data bytes
 * @param imageWidth
 *    Width of image in pixels
 * @param imageHeight
 *    Height of image in pixels
 * @param imageSlices
 *    Number of slices in image (1 if 2D)
 * @param dimensionType
 *    Texture dimension (2D or 3D)
 * @param pixelFormatType
 *    Pixel color components type
 * @param pixelOrigin
 *    Origin of first pixel in image
 * @param wrappingType
 *    Texture wrapping type
 * @param mipMappingType
 *    Mip mapping type
 * @param compressionType
 *    Type of compression
 * @param magnificationFilter
 *    Type of texture magnification
 * @param minificationFilter
 *    Type of texture minification
 * @param borderColor
 *    Color for texture border (used where there is no texture data)
 */
GraphicsTextureSettings::GraphicsTextureSettings(const uint8_t*        imageBytesPointer,
                                                 const int64_t         imageWidth,
                                                 const int64_t         imageHeight,
                                                 const int64_t         imageSlices,
                                                 const DimensionType   dimensionType,
                                                 const PixelFormatType pixelFormatType,
                                                 const PixelOrigin     pixelOrigin,
                                                 const WrappingType    wrappingType,
                                                 const MipMappingType  mipMappingType,
                                                 const CompressionType compressionType,
                                                 const GraphicsTextureMagnificationFilterEnum::Enum magnificationFilter,
                                                 const GraphicsTextureMinificationFilterEnum::Enum minificationFilter,
                                                 const std::array<float, 4>& borderColor)
: CaretObject(),
m_imageDataType(ImageDataType::POINTER),
m_imageBytesPointer(const_cast<uint8_t*>(imageBytesPointer)),
m_imageWidth(imageWidth),
m_imageHeight(imageHeight),
m_imageSlices(imageSlices),
m_dimensionType(dimensionType),
m_pixelFormatType(pixelFormatType),
m_pixelOrigin(pixelOrigin),
m_wrappingType(wrappingType),
m_mipMappingType(mipMappingType),
m_compressionType(compressionType),
m_magnificationFilter(magnificationFilter),
m_minificationFilter(minificationFilter),
m_borderColor(borderColor)
{
    switch (m_mipMappingType) {
        case MipMappingType::DISABLED:
        {
            AString filterName;
            switch (m_minificationFilter) {
                case GraphicsTextureMinificationFilterEnum::LINEAR:
                    break;
                case GraphicsTextureMinificationFilterEnum::LINEAR_MIPMAP_LINEAR:
                    filterName = "LINEAR_MIPMAP_LINEAR";
                    break;
                case GraphicsTextureMinificationFilterEnum::LINEAR_MIPMAP_NEAREST:
                    filterName = "LINEAR_MIPMAP_NEAREST";
                    break;
                case GraphicsTextureMinificationFilterEnum::NEAREST:
                    break;
                case GraphicsTextureMinificationFilterEnum::NEAREST_MIPMAP_LINEAR:
                    filterName = "NEAREST_MIPMAP_LINEAR";
                    break;
                case GraphicsTextureMinificationFilterEnum::NEAREST_MIPMAP_NEAREST:
                    filterName = "NEAREST_MIPMAP_NEAREST";
                    break;
            }
            
            if ( ! filterName.isEmpty()) {
                const AString msg("Mip mapping is DISABLED but the minification filter is set to a "
                                  "mip mapping type "
                                  + filterName
                                  + ".  This may result in nothing or just white displayed.  Either "
                                  "enable mip mapping or use a non-mip mapping minification filter.");
                CaretLogSevere(msg);
            }
        }
            break;
        case MipMappingType::ENABLED:
            break;
    }
}

/**
 * Destructor.
 */
GraphicsTextureSettings::~GraphicsTextureSettings()
{
}

/**
 * Copy constructor.
 * @param obj
 *    Object that is copied.
 */
GraphicsTextureSettings::GraphicsTextureSettings(const GraphicsTextureSettings& obj)
: CaretObject(obj)
{
    this->copyHelperGraphicsTextureSettings(obj);
}

/**
 * Assignment operator.
 * @param obj
 *    Data copied from obj to this.
 * @return 
 *    Reference to this object.
 */
GraphicsTextureSettings&
GraphicsTextureSettings::operator=(const GraphicsTextureSettings& obj)
{
    if (this != &obj) {
        CaretObject::operator=(obj);
        this->copyHelperGraphicsTextureSettings(obj);
    }
    return *this;    
}

/**
 * Helps with copying an object of this type.
 * @param obj
 *    Object that is copied.
 */
void 
GraphicsTextureSettings::copyHelperGraphicsTextureSettings(const GraphicsTextureSettings& obj)
{
    m_imageDataType       = obj.m_imageDataType;
    m_imageRgbaData       = obj.m_imageRgbaData;
    m_imageBytesPointer   = obj.m_imageBytesPointer;
    m_imageWidth          = obj.m_imageWidth;
    m_imageHeight         = obj.m_imageHeight;
    m_imageSlices         = obj.m_imageSlices;
    m_dimensionType       = obj.m_dimensionType;
    m_pixelFormatType     = obj.m_pixelFormatType;
    m_pixelOrigin         = obj.m_pixelOrigin;
    m_wrappingType        = obj.m_wrappingType;
    m_mipMappingType      = obj.m_mipMappingType;
    m_unpackAlignment     = obj.m_unpackAlignment;
    m_compressionType     = obj.m_compressionType;
    m_magnificationFilter = obj.m_magnificationFilter;
    m_minificationFilter  = obj.m_minificationFilter;
    m_borderColor         = obj.m_borderColor;
}

/**
 * @return Pointer to the qImage bytes
 */
const uint8_t*
GraphicsTextureSettings::getImageBytesPointer() const
{
    switch (m_imageDataType) {
        case ImageDataType::INVALID:
            CaretAssert(0);
            return NULL;
            break;
        case ImageDataType::POINTER:
            CaretAssert(m_imageBytesPointer);
            return m_imageBytesPointer;
            break;
        case ImageDataType::SHARED_PTR:
            CaretAssert(m_imageRgbaData.get());
            return m_imageRgbaData.get();
            break;
    }
    
    CaretAssert(0);
    return NULL;
}

/**
 * @return The widget of the image
 */
int64_t
GraphicsTextureSettings::getImageWidth() const
{
    return m_imageWidth;
}

/**
 * @return The height of the image
 */
int64_t
GraphicsTextureSettings::getImageHeight() const
{
    return m_imageHeight;
}

/**
 * @return The number of slices in the image (2D is 1)
 */
int64_t
GraphicsTextureSettings::getImageSlices() const
{
    return m_imageSlices;
}

/**
 * @return The dimension type
 */
GraphicsTextureSettings::DimensionType
GraphicsTextureSettings::getDimensionType() const
{
    return m_dimensionType;
}

/**
 * @return The pixel format
 */
GraphicsTextureSettings::PixelFormatType
GraphicsTextureSettings::getPixelFormatType() const
{
    return m_pixelFormatType;
}

/**
 * @return The pixel origin
 */
GraphicsTextureSettings::PixelOrigin
GraphicsTextureSettings::getPixelOrigin() const
{
    return m_pixelOrigin;
}

/**
 * @return The wrapping type
 */
GraphicsTextureSettings::WrappingType
GraphicsTextureSettings::getWrappingType() const
{
    return m_wrappingType;
}

/**
 * @return The mip mapping type
 */
GraphicsTextureSettings::MipMappingType
GraphicsTextureSettings::getMipMappingType() const
{
    return m_mipMappingType;
}

/**
 * @return The type of compression
 */
GraphicsTextureSettings::CompressionType
GraphicsTextureSettings::getCompressionType() const
{
    return m_compressionType;
}

/**
 * @return The magnification filter
 */
GraphicsTextureMagnificationFilterEnum::Enum
GraphicsTextureSettings::getMagnificationFilter() const
{
    return m_magnificationFilter;
}

/**
 * @return The unpack alignment.  See OpenGL documentation.  If wrong
 * image will be distorted (skewed).  May be needed if length of row in bytes
 * is not divisible by 4.
 */
int32_t
GraphicsTextureSettings::getUnpackAlignment() const
{
    return m_unpackAlignment;
}

void
GraphicsTextureSettings::setUnpackAlignmnet(const int32_t unpackAlignment)
{
    m_unpackAlignment = unpackAlignment;
}

/**
 * Set the magification filter
 * @param magificationFilter
 *     New value
 */
void
GraphicsTextureSettings::setMagnificationFilter(const GraphicsTextureMagnificationFilterEnum::Enum magnificationFilter)
{
    m_magnificationFilter = magnificationFilter;
}

/**
 * @return The minification filter
 */
GraphicsTextureMinificationFilterEnum::Enum
GraphicsTextureSettings::getMinificationFilter() const
{
    return m_minificationFilter;
}

/**
 * Set the minification filter
 * @param minificationFilter
 *     New value
 */
void
GraphicsTextureSettings::setMinificationFilter(const GraphicsTextureMinificationFilterEnum::Enum minificationFilter)
{
    m_minificationFilter = minificationFilter;
}

/**
 * @return The border color
 */
std::array<float, 4>
GraphicsTextureSettings::getBorderColor() const
{
    return m_borderColor;
}
/**
 * Get a description of this object's content.
 * @return String describing this object's content.
 */
AString 
GraphicsTextureSettings::toString() const
{
    AString txt;
    
    AString imagePointerTypeString("   Image Pointer Type: ");
    switch (m_imageDataType) {
        case ImageDataType::INVALID:
            imagePointerTypeString.append("INVALID");
            break;
        case ImageDataType::POINTER:
            imagePointerTypeString.append("POINTER");
            break;
        case ImageDataType::SHARED_PTR:
            imagePointerTypeString.append("SHARED_PTR");
            break;
    }
    txt.appendWithNewLine(imagePointerTypeString);
    
    AString dimensionTypeString("   Dimension Type: ");
    switch (m_dimensionType) {
        case DimensionType::NONE:
            dimensionTypeString.append("NONE");
            break;
        case DimensionType::FLOAT_STR_2D:
            dimensionTypeString.append("FLOAT_STR_2D");
            break;
        case DimensionType::FLOAT_STR_3D:
            dimensionTypeString.append("FLOAT_STR_3D");
            break;
    }
    txt.appendWithNewLine(dimensionTypeString);
    
    txt.appendWithNewLine("   Texture Width: "
                          + AString::number(m_imageWidth));
    txt.appendWithNewLine("   Texture Height: "
                          + AString::number(m_imageHeight));
    txt.appendWithNewLine("   Texture Slices: "
                          + AString::number(m_imageSlices));

    AString pixelFormatString("   Pixel Format Type: ");
    switch (m_pixelFormatType) {
        case PixelFormatType::NONE:
            pixelFormatString.append("NONE");
            break;
        case PixelFormatType::BGR:
            pixelFormatString.append("BGR");
            break;
        case PixelFormatType::BGRA:
            pixelFormatString.append("BGRA");
            break;
        case PixelFormatType::BGRX:
            pixelFormatString.append("BGRX");
            break;
        case PixelFormatType::RGB:
            pixelFormatString.append("RGB");
            break;
        case PixelFormatType::RGBA:
            pixelFormatString.append("RGBA");
            break;
    }
    txt.appendWithNewLine(pixelFormatString);
    
    AString pixelOriginString("   Pixel Origin: ");
    switch (m_pixelOrigin) {
        case PixelOrigin::NONE:
            pixelOriginString.append("NONE");
            break;
        case PixelOrigin::BOTTOM_LEFT:
            pixelOriginString.append("BOTTOM_LEFT");
            break;
        case PixelOrigin::TOP_LEFT:
            pixelOriginString.append("TOP_LEFT");
            break;
    }
    txt.appendWithNewLine(pixelOriginString);
    
    AString wrappingTypeString("   Wrapping Type: ");
    switch (m_wrappingType) {
        case WrappingType::CLAMP:
            wrappingTypeString.append("CLAMP");
            break;
        case WrappingType::CLAMP_TO_BORDER:
            wrappingTypeString.append("CLAMP_TO_BORDER");
            break;
        case WrappingType::REPEAT:
            wrappingTypeString.append("REPEAT");
            break;
    }
    txt.appendWithNewLine(wrappingTypeString);
    
    AString mipMappingTypeString("   Mip Mapping Type: ");
    switch (m_mipMappingType) {
        case MipMappingType::DISABLED:
            mipMappingTypeString.append("DISABLED");
            break;
        case MipMappingType::ENABLED:
            mipMappingTypeString.append("ENABLED");
            break;
    }
    txt.appendWithNewLine(mipMappingTypeString);
    
    AString compressionTypeString("   Compression Type: ");
    switch (m_compressionType) {
        case CompressionType::ENABLED:
            compressionTypeString.append("ENABLED");
            break;
        case CompressionType::DISABLED:
            compressionTypeString.append("DISABLED");
            break;
    }
    txt.appendWithNewLine(compressionTypeString);
    
    AString magnificationFilterString("   Magnification Filter: ");
    switch (m_magnificationFilter) {
        case GraphicsTextureMagnificationFilterEnum::LINEAR:
            magnificationFilterString.append("LINEAR");
            break;
        case GraphicsTextureMagnificationFilterEnum::NEAREST:
            magnificationFilterString.append("NEAREST");
            break;
    }
    txt.appendWithNewLine(magnificationFilterString);
    
    AString minificationFilterString("   Minification Filter: ");
    switch (m_minificationFilter) {
        case GraphicsTextureMinificationFilterEnum::NEAREST:
            minificationFilterString.append("NEAREST");
            break;
        case GraphicsTextureMinificationFilterEnum::LINEAR:
            minificationFilterString.append("LINEAR");
            break;
        case GraphicsTextureMinificationFilterEnum::LINEAR_MIPMAP_LINEAR:
            minificationFilterString.append("LINEAR_MIPMAP_LINEAR");
            break;
        case GraphicsTextureMinificationFilterEnum::LINEAR_MIPMAP_NEAREST:
            minificationFilterString.append("LINEAR_MIPMAP_NEAREST");
            break;
        case GraphicsTextureMinificationFilterEnum::NEAREST_MIPMAP_LINEAR:
            minificationFilterString.append("NEAREST_MIPMAP_LINEAR");
            break;
        case GraphicsTextureMinificationFilterEnum::NEAREST_MIPMAP_NEAREST:
            minificationFilterString.append("NEAREST_MIPMAP_NEAREST");
            break;
    }
    txt.appendWithNewLine(minificationFilterString);
    
    return txt;
}

