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

#include <cmath>
#include <cstdint>

#include <QBuffer>
#include <QColor>
#include <QDir>
#include <QImage>
#include <QImageReader>
#include <QImageWriter>
#include <QTime>

#define __IMAGE_FILE_DECLARE__
#include "ImageFile.h"
#undef __IMAGE_FILE_DECLARE__

#include "ApplicationInformation.h"
#include "BoundingBox.h"
#include "CaretAssert.h"
#include "CaretLogger.h"
#include "CaretPreferences.h"
#include "ControlPointFile.h"
#include "ControlPoint3D.h"
#include "DataFileException.h"
#include "DataFileContentInformation.h"
#include "DeveloperFlagsEnum.h"
#include "EventCaretPreferencesGet.h"
#include "EventManager.h"
#include "FileInformation.h"
#include "GiftiMetaData.h"
#include "GraphicsUtilitiesOpenGL.h"
#include "GraphicsPrimitiveV3fT2f.h"
#include "ImageCaptureDialogSettings.h"
#include "Matrix4x4.h"
#include "MathFunctions.h"
#include "RectangleTransform.h"
#include "SceneClass.h"
#include "UnitsConversion.h"
#include "VolumeFile.h"
#include "VolumeSpace.h"

using namespace caret;

static bool imageDebugFlag = false;

/**
 * Constructor.
 */
ImageFile::ImageFile()
: MediaFile(DataFileTypeEnum::IMAGE)
{
    initializeMembersImageFile();

    readFileMetaDataFromQImage();
}

/**
 * Initialize members of media file
 */
void
ImageFile::initializeMembersImageFile()
{
    m_controlPointFile.grabNew(new ControlPointFile());
    m_fileMetaData.grabNew(new GiftiMetaData());
    m_image = new QImage();
    m_graphicsPrimitive.reset();
    m_featuresImageGraphicsPrimitive.reset();
    m_pixelPrimitiveVertexStartIndex = -1;
    m_pixelPrimitiveVertexCount      = -1;
    m_planePrimitiveVertexStartIndex = -1;
    m_planePrimitiveVertexCount      = -1;
}

/**
 * Copy constructor.
 * @param imageFile
 *   Image file that is copied.
 */
ImageFile::ImageFile(const ImageFile& imageFile)
: MediaFile(imageFile)
{
    initializeMembersImageFile();
    
    if (m_image != NULL) {
        delete m_image;
    }
    if (imageFile.m_image != NULL) {
        m_image = new QImage(*imageFile.m_image);
    }
    else {
        m_image = new QImage();
    }
    
    m_fileMetaData.grabNew(new GiftiMetaData(*imageFile.m_fileMetaData));
    
    if (imageFile.m_controlPointFile) {
        /*CaretAssertMessage(0, "Need to implement copy contructor for ControlPointFile.");*/
        //m_controlPointFile.grabNew(new ControlPointFile(*imageFile.m_controlPointFile));
    }
    m_graphicsPrimitive.reset();
    m_featuresImageGraphicsPrimitive.reset();
    
    m_pixelPrimitiveVertexStartIndex = imageFile.m_pixelPrimitiveVertexStartIndex;;
    m_pixelPrimitiveVertexCount      = imageFile.m_pixelPrimitiveVertexCount;
    m_planePrimitiveVertexStartIndex = imageFile.m_planePrimitiveVertexStartIndex;
    m_planePrimitiveVertexCount      = imageFile.m_planePrimitiveVertexCount;
}

/**
 * Constructor that makes copy of QImage instance
 * @param qimage
 *    QImage that is copied to this image file.
 */
ImageFile::ImageFile(const QImage& qimage)
: MediaFile(DataFileTypeEnum::IMAGE)
{
    initializeMembersImageFile();
    if (m_image != NULL) {
        delete m_image;
    }
    m_image = new QImage(qimage);

    readFileMetaDataFromQImage();
}

/**
 * Constructs an image file from image data.
 *
 * @param imageDataRGBA
 *     Image data unsigned bytes with one byte for each
 *     red, green, blue, alpha.
 * @param imageWidth
 *     Width of image.
 * @param imageHeight
 *     Height of image.
 * @param imageOrigin
 *     Location of first pixel in the image data.
 */
ImageFile::ImageFile(const unsigned char* imageDataRGBA,
                     const int imageWidth,
                     const int imageHeight,
                     const IMAGE_DATA_ORIGIN_LOCATION imageOrigin)
: MediaFile(DataFileTypeEnum::IMAGE)
{
    initializeMembersImageFile();

    if (m_image != NULL) {
        delete m_image;
    }
    m_image = new QImage(imageWidth,
                             imageHeight,
                             QImage::Format_ARGB32);
    readFileMetaDataFromQImage();
    
    bool isOriginAtTop = false;
    switch (imageOrigin) {
        case IMAGE_DATA_ORIGIN_AT_BOTTOM:
            isOriginAtTop = false;
            break;
        case IMAGE_DATA_ORIGIN_AT_TOP:
            isOriginAtTop = true;
            break;
    }
    
    /*
     * Documentation for QImage states that setPixel may be very costly
     * and recommends using the scanLine() method to access pixel data.
     */
    for (int y = 0; y < imageHeight; y++) {
        const int scanLineIndex = (isOriginAtTop
                                   ? y
                                   : imageHeight -y - 1);
        QRgb* rgbScanLine = (QRgb*)m_image->scanLine(scanLineIndex);
        
        for (int x = 0; x < imageWidth; x++) {
            const int32_t contentOffset = (((y * imageWidth) * 4)
                                           + (x * 4));
            const int red   = imageDataRGBA[contentOffset];
            const int green = imageDataRGBA[contentOffset+1];
            const int blue  = imageDataRGBA[contentOffset+2];
            const int alpha  = imageDataRGBA[contentOffset+3];
            QColor color(red,
                         green,
                         blue,
                         alpha);
            
            QRgb* pixel = &rgbScanLine[x];
            *pixel = color.rgba();
        }
    }
    readFileMetaDataFromQImage();
}

/**
 * Constructs an image file with black and zero alpha pixels
 *
 * @param imageWidth
 *     Width of image.
 * @param imageHeight
 *     Height of image.
 */
ImageFile::ImageFile(const int imageWidth,
                     const int imageHeight)
: MediaFile(DataFileTypeEnum::IMAGE)
{
    initializeMembersImageFile();
    
    if (m_image != NULL) {
        delete m_image;
    }
    m_image = new QImage(imageWidth,
                         imageHeight,
                         QImage::Format_ARGB32);
    readFileMetaDataFromQImage();
    
    /*
     * fill image with black with zero alpha
     */
    m_image->fill(QColor(0, 0, 0, 0));

    readFileMetaDataFromQImage();
}

/**
 * Destructor.
 */
ImageFile::~ImageFile()
{
    if (m_image != NULL) {
        delete m_image;
        m_image = NULL;
    }
}

/**
 * Clears current file data in memory.
 */
void
ImageFile::clear()
{
    if (m_image != NULL) {
        delete m_image;
    }
    initializeMembersImageFile();

    readFileMetaDataFromQImage();
    this->clearModified();
}

/**
 * Clear this file's modified status
 */
void
ImageFile::clearModified()
{
    MediaFile::clearModified();
    m_fileMetaData->clearModified();
}

/**
 * @return True if this file is modified, else falsel
 */
bool
ImageFile::isModified() const
{
    if (MediaFile::isModified()) {
        return true;
    }
    if (m_fileMetaData->isModified()) {
        return true;
    }
    return false;
}

/**
 * @return Number of frames in the file
 */
int32_t
ImageFile::getNumberOfFrames() const
{
    if (m_image != NULL) {
        return 1;
    }
    return 0;
}

/**
 * @return Name of frame at given index.
 * @param frameIndex Index of the frame
 */
AString
ImageFile::getFrameName(const int32_t frameIndex) const
{
    AString frameName;
    if (frameName.isEmpty()) {
        frameName = MediaFile::getFrameName(frameIndex);
    }
    
    return frameName;
}

/**
 * @return Get access to the file's metadata.
 */
GiftiMetaData*
ImageFile::getFileMetaData()
{
    return m_fileMetaData;
}

/**
 * @return Get access to unmodifiable file's metadata.
 */
const GiftiMetaData*
ImageFile::getFileMetaData() const
{
    return m_fileMetaData;
}

/**
 * @return true if the file is is empty (image contains no pixels).
 */

bool
ImageFile::isEmpty() const
{
    return (m_image->width() <= 0);
}

/**
 * @return A pointer to the QImage in this file.
 */
const QImage*
ImageFile::getAsQImage() const
{
    writeFileMetaDataToQImage();
    return m_image;
}

/**
 * Set the image in this file from a QImage.
 * @param qimage
 *    Image that is copied to this file.
 */
void
ImageFile::setFromQImage(const QImage& qimage)
{
    if (m_image != NULL) {
        delete m_image;
    }
    m_image = new QImage(qimage);
    readFileMetaDataFromQImage();
    this->setModified();
}

/**
 * Set the dots per meter.
 *
 * @param x
 *    Dots per meter for X dimension.
 * @param y
 *    Dots per meter for Y dimension.
 */
void
ImageFile::setDotsPerMeter(const int x,
                           const int y)
{
    m_image->setDotsPerMeterX(x);
    m_image->setDotsPerMeterY(y);
}


/**
 * Examines the image to find the rectangular of the object in the image
 * by examining pixels in the background color.
 * @param backgroundColor
 *    RGB components range 0-255.
 * @param objectBoundsOut
 *    4-dimensional array containing the region that excludes
 *    the backround around the image's object.
 */
void
ImageFile::findImageObject(const uint8_t backgroundColor[3],
                           int objectBoundsOut[4]) const
{
    //
    // Dimensions of image
    //
    const int numX = m_image->width();
    const int numY = m_image->height();
    
    //
    // Initialize output
    //
    objectBoundsOut[0] = 0;
    objectBoundsOut[1] = 0;
    objectBoundsOut[2] = numX - 1;
    objectBoundsOut[3] = numY - 1;
    
    //
    // Find left
    //
    bool gotPixelFlag = false;
    for (int i = 0; i < numX; i++) {
        for (int j = 0; j < numY; j++) {
            const QRgb pixel = m_image->pixel(i, j);
            if ((qRed(pixel) != backgroundColor[0]) ||
                (qGreen(pixel) != backgroundColor[1]) ||
                (qBlue(pixel)  != backgroundColor[2])) {
                objectBoundsOut[0] = i;
                gotPixelFlag = true;
                break;
            }
        }
        if (gotPixelFlag) {
            break;
        }
    }
    
    //
    // Find right
    //
    gotPixelFlag = false;
    for (int i = (numX - 1); i >= 0; i--) {
        for (int j = 0; j < numY; j++) {
            const QRgb pixel = m_image->pixel(i, j);
            if ((qRed(pixel) != backgroundColor[0]) ||
                (qGreen(pixel) != backgroundColor[1]) ||
                (qBlue(pixel)  != backgroundColor[2])) {
                objectBoundsOut[2] = i;
                gotPixelFlag = true;
                break;
            }
        }
        if (gotPixelFlag) {
            break;
        }
    }
    
    //
    // Find top
    //
    gotPixelFlag = false;
    for (int j = 0; j < numY; j++) {
        for (int i = 0; i < numX; i++) {
            const QRgb pixel = m_image->pixel(i, j);
            if ((qRed(pixel) != backgroundColor[0]) ||
                (qGreen(pixel) != backgroundColor[1]) ||
                (qBlue(pixel)  != backgroundColor[2])) {
                objectBoundsOut[1] = j;
                gotPixelFlag = true;
                break;
            }
        }
        if (gotPixelFlag) {
            break;
        }
    }
    
    //
    // Find bottom
    //
    gotPixelFlag = false;
    for (int j = (numY - 1); j >= 0; j--) {
        for (int i = 0; i < numX; i++) {
            const QRgb pixel = m_image->pixel(i, j);
            if ((qRed(pixel) != backgroundColor[0]) ||
                (qGreen(pixel) != backgroundColor[1]) ||
                (qBlue(pixel)  != backgroundColor[2])) {
                objectBoundsOut[3] = j;
                gotPixelFlag = true;
                break;
            }
        }
        if (gotPixelFlag) {
            break;
        }
    }
}

/**
 * Add a margin to this image.
 * @param marginSize
 *    Number of pixels in the margin.
 * @param backgroundColor used for the added pixels.
 *    RGB components range 0-255.
 */
void
ImageFile::addMargin(const int marginSize,
                     const uint8_t backgroundColor[3])
{
    this->addMargin(marginSize, marginSize, backgroundColor);
}

/**
 * Add a margin to this image.
 * @param image
 *    Image to which margin is added.
 * @param marginSizeX
 *    Number of pixels in the margin along x-axis.
 * @param marginSizeY
 *    Number of pixels in the margin along y-axis.
 * @param backgroundColor used for the added pixels.
 *    RGB components range 0-255.
 */
void
ImageFile::addMargin(const int marginSizeX,
                     const int marginSizeY,
                     const uint8_t backgroundColor[3])
{
    if ((marginSizeX <= 0) && (marginSizeY <= 0)) {
        return;
    }
    
    //
    // Add margin
    //
    const int width = m_image->width();
    const int height = m_image->height();
    const int newWidth = width + marginSizeX * 2;
    const int newHeight = height + marginSizeY * 2;
    QRgb backgroundColorRGB = qRgba(backgroundColor[0],
                                    backgroundColor[1],
                                    backgroundColor[2],
                                    255);
    
    GiftiMetaData fileMetaDataCopy(*m_fileMetaData);
    
    //
    // Insert image
    //
    ImageFile imageFile;
    imageFile.setFromQImage(QImage(newWidth, newHeight, m_image->format()));
    imageFile.m_image->fill(backgroundColorRGB);
    try {
        imageFile.insertImage(*m_image, marginSizeX, marginSizeY);
        this->setFromQImage(*imageFile.getAsQImage());
        
        /*
         * Preserve metadata
         */
        *m_fileMetaData = fileMetaDataCopy;
    }
    catch (DataFileException& e) {
        CaretLogWarning(e.whatString());
    }
    this->setModified();
}

/**
 * Crop an image by removing the background from the object in the image
 * but keeping a margin of the given size around the image.
 * @param marginSize
 *    Number of pixels in the margin around the image's object.
 * @param backgroundColor
 *    Color of background that ranges 0-255.
 */
void
ImageFile::cropImageRemoveBackground(const int marginSize,
                                     const uint8_t backgroundColor[3])
{
    //
    // Get cropping bounds
    //
    int leftTopRightBottom[4];
    this->findImageObject(backgroundColor,
                          leftTopRightBottom);
    CaretLogFine("cropping: "
                   + AString::fromNumbers(leftTopRightBottom, 4, " "));
    
    const int currentWidth = m_image->width();
    const int currentHeight = m_image->height();
    
    //
    // If cropping is valid
    //
    const int width = leftTopRightBottom[2] - leftTopRightBottom[0] + 1;
    const int height = leftTopRightBottom[3] - leftTopRightBottom[1] + 1;
    if ((width != currentWidth)
        || (height != currentHeight)) {
        if ((width > 1) &&
            (height > 1)) {
            QImage copyImage = this->getAsQImage()->copy(leftTopRightBottom[0],
                                                         leftTopRightBottom[1],
                                                         width,
                                                         height);
            if ( ! copyImage.isNull()) {
                if ((copyImage.width() > 0)
                    && (copyImage.height() > 0)) {
                    GiftiMetaData fileMetaDataCopy(*m_fileMetaData);
                    this->setFromQImage(copyImage);
                    *m_fileMetaData = fileMetaDataCopy;
                }
            }
            
            this->setModified();
        }
    }
    
    //
    // Process margin
    //
    if (marginSize > 0) {
        this->addMargin(marginSize,
                        backgroundColor);
    }
}

/**
 * Replace the contents of this image by combining the input images and
 * retaining aspect and stretching and filling if needed.
 * @param imageFiles
 *    Images that are combined.
 * @param numImagesPerRow
 *    Number of images in each row.
 * @param backgroundColor
 *    Color of background that ranges 0-255.
 */
void
ImageFile::combinePreservingAspectAndFillIfNeeded(const std::vector<ImageFile*>& imageFiles,
                                                  const int numImagesPerRow,
                                                  const uint8_t backgroundColor[3])
{
    const int numImages = static_cast<int>(imageFiles.size());
    if (numImages <= 0) {
        return;
    }
    if (numImages == 1) {
        this->setFromQImage(*imageFiles[0]->m_image);
        return;
    }
    
    QRgb backgroundColorRGB = qRgba(backgroundColor[0],
                                    backgroundColor[1],
                                    backgroundColor[2],
                                    0);
    //
    // Resize all images but do not stretch
    // need to retain aspect ratio but all must
    // be the same size in X & Y
    //
    
    //
    // Find max width and height of input images
    //
    int maxImageWidth = 0;
    int maxImageHeight = 0;
    for (int i = 0; i < numImages; i++) {
        //
        // Track max width/height
        //
        maxImageWidth = std::max(maxImageWidth, imageFiles[i]->m_image->width());
        maxImageHeight = std::max(maxImageHeight, imageFiles[i]->m_image->height());
    }
    
    //
    // Compute size of output image and create it
    //
    const int outputImageSizeX = maxImageWidth * numImagesPerRow;
    const int numberOfRows = (numImages / numImagesPerRow)
    + (((numImages % numImagesPerRow) != 0) ? 1 : 0);
    const int outputImageSizeY = maxImageHeight * numberOfRows;
    QImage combinedImage(outputImageSizeX,
                         outputImageSizeY,
                         imageFiles[0]->m_image->format());
    combinedImage.fill(backgroundColorRGB);
    
    
    //
    // Loop through the images
    //
    int rowCounter = 0;
    int columnCounter = 0;
    for (int i = 0; i < numImages; i++) {
        //
        // Scale image
        //
        const QImage imageScaled = imageFiles[i]->m_image->scaled(maxImageWidth,
                                                                maxImageHeight,
                                                                Qt::KeepAspectRatio,
                                                                Qt::SmoothTransformation);
        
        //
        // Compute position of where image should be inserted
        //
        const int marginX = (maxImageWidth - imageScaled.width()) / 2;
        const int marginY = (maxImageHeight - imageScaled.height()) / 2;
        const int positionX = columnCounter * maxImageWidth + marginX;
        const int positionY = rowCounter * maxImageHeight + marginY;
        
        //
        // Insert into output image
        //
        try {
            ImageFile::insertImage(imageScaled,
                                   combinedImage,
                                   positionX,
                                   positionY);
        }
        catch (DataFileException& e) {
            CaretLogWarning("QImageFile::insertImage() error: "
                            + e.whatString());
        }
        
        //
        // Update row and column counters
        //
        columnCounter++;
        if (columnCounter >= numImagesPerRow) {
            columnCounter = 0;
            rowCounter++;
        }
    }
    
    this->setFromQImage(combinedImage);
}

/**
 * Read the image file.
 * @param filename
 *    Name of image file.
 * @throws DataFileException
 *    If error reading image.
 */
void
ImageFile::readFile(const AString& filename)
{
    clear();
    
    checkFileReadability(filename);
    
    if (filename.isEmpty()) {
        throw DataFileException(filename + "Filename for reading is isEmpty");
    }
    
    this->setFileName(filename);
    
    if ( ! m_image->load(filename)) {
        clear();
        throw DataFileException(filename + "Unable to load file.");
    }
    
    m_image = limitImageDimensions(m_image,
                                   filename);
    
    /*
     * Format must be RGB or ARGB for compatibility with OpenGL
     */
    verifyFormatCompatibleWithOpenGL();

    readFileMetaDataFromQImage();
    
    this->clearModified();
}

/**
 * Limit the dimensions of an the image
 * @param image
 *    Image that is limited in dimensions
 * @param filename
 *    Name of image file
 * @return
 *    Original image or new image that is resized version of original image
 */
QImage*
ImageFile::limitImageDimensions(QImage* image,
                                const AString& filename)
{
    QImage* imageOut(image);
    
    if ( ! image->isNull()) {
        switch (ApplicationInformation::getApplicationType()) {
            case ApplicationTypeEnum::APPLICATION_TYPE_COMMAND_LINE:
                break;
            case ApplicationTypeEnum::APPLICATION_TYPE_GRAPHICAL_USER_INTERFACE:
            {
                const int32_t width(image->width());
                const int32_t height(image->height());
                const int32_t maxDim(GraphicsUtilitiesOpenGL::getTextureWidthHeightMaximumDimension());
                if ((width > maxDim)
                    || (height > maxDim)) {
                    QImage newImage;
                    if (width > height) {
                        newImage = image->scaledToWidth(maxDim);
                    }
                    else {
                        newImage = image->scaledToHeight(maxDim);
                    }
                    if ( ! newImage.isNull()) {
                        delete image;
                        imageOut = new QImage(newImage);
                        CaretLogWarning("Rescaled image "
                                        + filename
                                        + " from size ("
                                        + AString::number(width)
                                        + ", "
                                        + AString::number(height)
                                        + ") to ("
                                        + AString::number(imageOut->width())
                                        + ", "
                                        + AString::number(imageOut->height())
                                        + ")");
                    }
                }
            }
                break;
            case ApplicationTypeEnum::APPLICATION_TYPE_INVALID:
                break;
        }
    }
    
    return imageOut;
}

/**
 * Insert an image into this image which must be large enough for insertion of image.
 * @param otherImage
 *    Image that is inserted into this image.
 * @param x
 *    X position of where image is inserted.
 * @param y
 *    Y position of where image is inserted.
 * @throws DataFileException
 *    If error inserting image.
 */
void
ImageFile::insertImage(const QImage& otherImage,
                       const int x,
                       const int y)
{
    ImageFile::insertImage(otherImage,
                           *m_image,
                           x,
                           y);
    this->setModified();
}

/**
 * insert an image into another image.
 *
 * Insert an image into another image which must be large enough for insertion of image.
 * @param insertThisImage
 *    Image that is inserted into other image.
 * @param intoThisImage
 *    Image that receives insertions of other image.
 * @param x
 *    X position of where image is inserted.
 * @param y
 *    Y position of where image is inserted.
 * @throws DataFileException
 *    If error inserting image.
 */
void
ImageFile::insertImage(const QImage& insertThisImage,
                       QImage& intoThisImage,
                       const int positionX,
                       const int positionY)
{
    if (positionX < 0) {
        throw DataFileException("X position is less than zero.");
    }
    if (positionY < 0) {
        throw DataFileException("Y position is less than zero.");
    }
    
    const int otherWidth = insertThisImage.width();
    const int otherHeight = insertThisImage.height();
    
    const int myWidth = intoThisImage.width();
    const int myHeight = intoThisImage.height();
    
    if ((otherWidth + positionX) > myWidth) {
        throw DataFileException("This image is not large enough to insert other image.");
    }
    if ((otherHeight + positionY) > myHeight) {
        throw DataFileException("This image is not large enough to insert other image.");
    }
    
    for (int i = 0; i < otherWidth; i++) {
        for (int j = 0; j < otherHeight; j++) {
            intoThisImage.setPixel(positionX + i,
                                   positionY + j,
                                   insertThisImage.pixel(i, j));
        }
    }
}

/**
 * Scale the given image to the given width and height while preserving
 * the aspect ratio.  If the
 *
 * @param image
 *     The image
 * @param width
 *     Width of the image
 * @param height
 *     Height of the image
 * @param fillColor
 *     If not NULL, padded region is this color
 * @return
 *     Image that will be the requested width and height or
 *     a null image (.isNull()) if error.
 */
QImage
ImageFile::scaleToSizeWithPadding(const QImage& image,
                                  const int width,
                                  const int height,
                                  const QColor* fillColor)
{
    /*
     * Invalid image tests
     */
    if (image.isNull()) {
        return image;
    }
    if ((image.width() <= 0)
        || (image.height() <= 0)) {
        return QImage();
    }
    
    /*
     * Nothing to do if image is correct size
     */
    if ((image.width() == width)
        && (image.height() == height)) {
        return image;
    }
    
    const QImage scaledImage = image.scaled(width,
                                      height,
                                      Qt::KeepAspectRatio,
                                      Qt::SmoothTransformation);
    const int scaledWidth  = scaledImage.width();
    const int scaledHeight = scaledImage.height();
    if ((scaledWidth == width)
        && (scaledHeight == height)) {
        return scaledImage;
    }
    else if (scaledWidth > width) {
        CaretLogSevere("Image scale width was made larger="
                       + QString::number(scaledWidth)
                       + " than requested="
                       + QString::number(width));
        return QImage();
    }
    else if (scaledHeight > height) {
        CaretLogSevere("Image scale height was made larger="
                       + QString::number(scaledHeight)
                       + " than requested="
                       + QString::number(height));
        return QImage();
    }

    QImage outputImage(width,
                       height,
                       image.format());
    if (fillColor != NULL) {
        outputImage.fill(*fillColor);
    }
    else {
        outputImage.fill(Qt::black);
    }
    
    const int insertX = (width - scaledWidth) / 2;
    const int insertY = (height - scaledHeight) / 2;
    
    try {
        ImageFile::insertImage(scaledImage,
                               outputImage,
                               insertX,
                               insertY);
    }
    catch (const DataFileException& dfe) {
        CaretLogSevere(dfe.whatString());
        outputImage = QImage();
    }
    
    return outputImage;
}

/**
 * Compare a file for unit testing (tolerance ignored).
 *
 * @param dataFile
 *    Data files that is compared to this data file.
 * @param tolerance
 *    Allowable difference at each pixel.
 * @param messageOut
 *    Message describing differences.
 * @return
 *    True if files are within tolerance, else false.
 */
bool
ImageFile::compareFileForUnitTesting(const DataFile* dataFile,
                                     const float tolerance,
                                     AString& messageOut) const
{
    //
    // Cast to an image file
    //
    const ImageFile* img = dynamic_cast<const ImageFile*>(dataFile);
    if (img == NULL) {
        messageOut = ("ERROR: File for comparison ("
                      + dataFile->getFileName()
                      + " does not appear to be an image file.");
        return false;
    }
    
    //
    // Get the image from the other file
    //
    const QImage* otherImage = img->getAsQImage();
    
    //
    // Confirm width/height
    //
    const int width = m_image->width();
    const int height = m_image->height();
    if ((width != otherImage->width()) ||
        (height != otherImage->height())) {
        messageOut = "The images are of different height and/or width.";
        return false;
    }
    
    //
    // compare pixels
    //
    int pixelCount = 0;
    for (int i = 0; i < width; i++) {
        for (int j = 0; j < height; j++) {
            QColor im1 = m_image->pixel(i, j);
            QColor im2 = otherImage->pixel(i, j);
            if ((abs(im1.red() - im2.red()) > tolerance) ||
                (abs(im1.green() - im2.green()) > tolerance) ||
                (abs(im1.blue() - im2.blue()) > tolerance)) {
                pixelCount++;
            }
        }
    }
    
    if (pixelCount > 0) {
        const float pct = static_cast<float>(pixelCount * 100.0)
        / static_cast<float>(width * height);
        messageOut = QString::number(pct, 'f', 2)
        + "% pixels in the image do not match.";
        return false;
    }
    
    return true;
}

/**
 * Write the image file.
 * @param filename
 *    Name of image file.
 * @throws DataFileException
 *    If error writing image.
 */
void
ImageFile::writeFile(const AString& filename)
{
    if ( ! DataFileTypeEnum::isValidWriteFileExtension(filename,
                                                       DataFileTypeEnum::IMAGE)) {
        throw DataFileException(filename,
                                "Filename's extension does not match a image file type that is writable.");
    }
    
    checkFileWritability(filename);
    
    this->setFileName(filename);
    
    AString errorMessage;
    if (m_image->width() <= 0) {
        errorMessage = "Image width is zero.";
    }
    if (m_image->height() <= 0) {
        if ( ! errorMessage.isEmpty()) errorMessage += "\n";
        errorMessage = "Image height is zero.";
    }
    if ( ! errorMessage.isEmpty()) {
        throw DataFileException(filename + "  " + errorMessage);
    }
    
    FileInformation fileInfo(this->getFileName());
    AString format = fileInfo.getFileExtension().toUpper();
    if (format == "JPG") {
        format = "JPEG";
    }
    
    QImageWriter writer(filename, format.toLatin1());
    if (writer.supportsOption(QImageIOHandler::Quality)) {
        if (format.compare("png", Qt::CaseInsensitive) == 0) {
            const int quality = 1;
            writer.setQuality(quality);
        }
        else {
            const int quality = 100;
            writer.setQuality(quality);
        }
    }
    
    if (writer.supportsOption(QImageIOHandler::CompressionRatio) && writer.compression() == 0) {
        writer.setCompression(1);
    }
    
    writeFileMetaDataToQImage();
    
    if ( ! writer.write(*m_image)) {
        throw DataFileException(writer.errorString());
    }
    
    this->clearModified();
}

/**
 * Get the image file filters for the supported image types for saving image files
 *
 * @param imageFileFilters
 *    Output filled with the filters for supported image types.
 * @param defaultFilter
 *    Filter for the preferred image type.
 */
void
ImageFile::getSaveQFileDialogImageFilters(std::vector<AString>& imageFileFilters,
                                          AString& defaultFilter)
{
    DataFileTypeEnum::getSaveQFileDialogImageFilters(imageFileFilters,
                                                     defaultFilter);
}

/**
 * Resize the image to the given width while preserving the aspect ratio
 * of the image.
 *
 * @param width
 *     Width for image.
 */
void
ImageFile::resizeToWidth(const int32_t width)
{
    GiftiMetaData fileMetaDataCopy(*m_fileMetaData);
    
    CaretAssert(m_image);
    *m_image = m_image->scaledToWidth(width,
                                      Qt::SmoothTransformation);
    
    *m_fileMetaData = fileMetaDataCopy;
}

/**
 * Resize the image to the given height while preserving the aspect ratio
 * of the image.
 *
 * @param height
 *     Height for image.
 */
void
ImageFile::resizeToHeight(const int32_t height)
{
    GiftiMetaData fileMetaDataCopy(*m_fileMetaData);
    
    CaretAssert(m_image);
    *m_image = m_image->scaledToHeight(height,
                                       Qt::SmoothTransformation);
    
    *m_fileMetaData = fileMetaDataCopy;
}

/**
 * Resize the image so that its width is no larger than the given value.
 * If the image's current width is less than the given value, no
 * resizing takes place.
 *
 * @param maximumWidth
 *     Maximum width for the image.
 */
void
ImageFile::resizeToMaximumWidth(const int32_t maximumWidth)
{
    CaretAssert(m_image);
    const int32_t width = m_image->width();
    
    if (width > maximumWidth) {
        GiftiMetaData fileMetaDataCopy(*m_fileMetaData);
        *m_image = m_image->scaledToWidth(maximumWidth,
                                          Qt::SmoothTransformation);
        *m_fileMetaData = fileMetaDataCopy;
    }
}

/**
 * Resize the image so that its height is no larger than the given value.
 * If the image's current height is less than the given value, no
 * resizing takes place.
 *
 * @param maximumHeight
 *     Maximum height for the image.
 */
void
ImageFile::resizeToMaximumHeight(const int32_t maximumHeight)
{
    CaretAssert(m_image);
    const int32_t height = m_image->height();
    
    if (height > maximumHeight) {
        GiftiMetaData fileMetaDataCopy(*m_fileMetaData);
        *m_image = m_image->scaledToHeight(maximumHeight,
                                          Qt::SmoothTransformation);
        *m_fileMetaData = fileMetaDataCopy;
    }
}

/**
 * Resize the image so that its maximum dimension is the given value
 * yet preserves the aspect ratio of the image.  If the maximum dimension
 * is less than the given value, no resizing takes place.
 *
 * @param maximumWidthOrHeight
 *     Maximum dimension for the image.
 */
void
ImageFile::resizeToMaximumWidthOrHeight(const int32_t maximumWidthOrHeight)
{
    CaretAssert(m_image);
    
    const int32_t width = m_image->width();
    const int32_t height = m_image->height();
    
    if ((width > 0)
        && (height > 0)) {
        if (width > height) {
            resizeToMaximumWidth(maximumWidthOrHeight);
        }
        else {
            if (height > maximumWidthOrHeight) {
                resizeToMaximumWidth(maximumWidthOrHeight);
            }
        }
    }
}

/**
 * Get the RGBA bytes from the image.
 *
 * @param bytesRGBA
 *    The RGBA bytes in the image.
 * @param widthOut
 *    Width of the image.
 * @param heightOut
 *    Height of the image.
 * @param imageOrigin
 *     Location of first pixel in the image data.
 * @return
 *    True if the bytes, width, and height are valid, else false.
 */
bool
ImageFile::getImageBytesRGBA(const IMAGE_DATA_ORIGIN_LOCATION imageOrigin,
                             std::vector<uint8_t>& bytesRGBA,
                             int32_t& widthOut,
                             int32_t& heightOut) const
{
    bytesRGBA.clear();
    widthOut = 0;
    heightOut = 0;
    
    if (m_image != NULL) {
        widthOut  = m_image->width();
        heightOut = m_image->height();
        if ((widthOut > 0)
            && (heightOut > 0)) {
            
            bytesRGBA.resize(widthOut * heightOut * 4);
            
            bool isOriginAtTop = false;
            switch (imageOrigin) {
                case IMAGE_DATA_ORIGIN_AT_BOTTOM:
                    isOriginAtTop = false;
                    break;
                case IMAGE_DATA_ORIGIN_AT_TOP:
                    isOriginAtTop = true;
                    break;
            }
            
            /*
             * Documentation for QImage states that setPixel may be very costly
             * and recommends using the scanLine() method to access pixel data.
             */
            for (int64_t y = 0; y < heightOut; y++) {
                const int64_t scanLineIndex = (isOriginAtTop
                                           ? y
                                           : heightOut -y - 1);
                const uchar* scanLine = m_image->scanLine(scanLineIndex);
                QRgb* rgbScanLine = (QRgb*)scanLine;
                
                for (int64_t x = 0; x < widthOut; x++) {
                    const int64_t contentOffset = (((y * widthOut) * 4)
                                                   + (x * 4));
                    QRgb& rgb = rgbScanLine[x];
                    bytesRGBA[contentOffset] = static_cast<uint8_t>(qRed(rgb));
                    bytesRGBA[contentOffset+1] = static_cast<uint8_t>(qGreen(rgb));
                    bytesRGBA[contentOffset+2] = static_cast<uint8_t>(qBlue(rgb));
                    bytesRGBA[contentOffset+3] = static_cast<uint8_t>(qAlpha(rgb));
                }
            }
            return true;
        }
    }
    
    return false;
}

/**
 * Get the RGBA bytes from the given QImage.
 *
 * @parm qImage
 *    The QImage
 * @param bytesRGBA
 *    The RGBA bytes in the image.
 * @param widthOut
 *    Width of the image.
 * @param heightOut
 *    Height of the image.
 * @param imageOrigin
 *     Location of first pixel in the image data.
 * @return
 *    True if the bytes, width, and height are valid, else false.
 */
bool
ImageFile::getImageBytesRGBA(const QImage* qImage,
                              const IMAGE_DATA_ORIGIN_LOCATION imageOrigin,
                              std::vector<uint8_t>& bytesRGBA,
                              int32_t& widthOut,
                              int32_t& heightOut)
{
    CaretAssert(qImage);
    bytesRGBA.clear();
    widthOut = 0;
    heightOut = 0;
    
    if (qImage != NULL) {
        widthOut  = qImage->width();
        heightOut = qImage->height();
        if ((widthOut > 0)
            && (heightOut > 0)) {
            
            bytesRGBA.resize(widthOut * heightOut * 4);
            
            bool isOriginAtTop = false;
            switch (imageOrigin) {
                case IMAGE_DATA_ORIGIN_AT_BOTTOM:
                    isOriginAtTop = false;
                    break;
                case IMAGE_DATA_ORIGIN_AT_TOP:
                    isOriginAtTop = true;
                    break;
            }
            
            /*
             * Documentation for QImage states that setPixel may be very costly
             * and recommends using the scanLine() method to access pixel data.
             */
            for (int64_t y = 0; y < heightOut; y++) {
                const int64_t scanLineIndex = (isOriginAtTop
                                               ? y
                                               : heightOut -y - 1);
                const uchar* scanLine = qImage->scanLine(scanLineIndex);
                QRgb* rgbScanLine = (QRgb*)scanLine;
                
                for (int64_t x = 0; x < widthOut; x++) {
                    const int64_t contentOffset = (((y * widthOut) * 4)
                                                   + (x * 4));
                    QRgb& rgb = rgbScanLine[x];
                    bytesRGBA[contentOffset] = static_cast<uint8_t>(qRed(rgb));
                    bytesRGBA[contentOffset+1] = static_cast<uint8_t>(qGreen(rgb));
                    bytesRGBA[contentOffset+2] = static_cast<uint8_t>(qBlue(rgb));
                    bytesRGBA[contentOffset+3] = static_cast<uint8_t>(qAlpha(rgb));
                }
            }
            return true;
        }
    }
    
    return false;
}

/**
 * Get the pixel RGBA at the given pixel I and J.
 *
 * @param tabIndex
 *    Index of the tab.
 * @param overlayIndex
 *    Index of overlay
 * @param pixelLogicalIndex
 *     Logical pixel index that may not have (0, 0) in top left corner
 * @param pixelRGBAOut
 *     RGBA at Pixel I, J
 * @return
 *     True if valid, else false.
 */
bool
ImageFile::getPixelRGBA(const int32_t /*tabIndex*/,
                        const int32_t /*overlayIndex*/,
                        const PixelLogicalIndex& pixelLogicalIndex,
                        uint8_t pixelRGBAOut[4]) const
{
    if (m_image != NULL) {
        const int32_t w = m_image->width();
        const int32_t h = m_image->height();
        
        const int64_t pixelI(pixelLogicalIndex.getI());
        const int64_t pixelJ(pixelLogicalIndex.getJ());
        if ((pixelI >= 0)
            && (pixelI < w)
            && (pixelJ >= 0)
            && (pixelJ < h)) {
            const QRgb rgb = m_image->pixel(pixelI,
                                            pixelJ);
            pixelRGBAOut[0] = static_cast<uint8_t>(qRed(rgb));
            pixelRGBAOut[1] = static_cast<uint8_t>(qGreen(rgb));
            pixelRGBAOut[2] = static_cast<uint8_t>(qBlue(rgb));
            pixelRGBAOut[3] = static_cast<uint8_t>(qAlpha(rgb));
            
            return true;
        }
    }
    
    return false;

}

/**
 * Set the pixel RGBA at the given pixel I and J.
 *
 * @param tabIndex
 *    Index of the tab.
 * @param overlayIndex
 *    Index of overlay
 * @param pixelLogicalIndex
 *     Logical pixel index that may not have (0, 0) in top left corner
 * @param pixelRGBAOut
 *     RGBA at Pixel I, J
 */
void
ImageFile::setPixelRGBA(const int32_t /*tabIndex*/,
                        const int32_t /*overlayIndex*/,
                        const PixelLogicalIndex& pixelLogicalIndex,
                        const uint8_t pixelRGBA[4])
{
    if (m_image != NULL) {
        const int64_t pixelI(pixelLogicalIndex.getI());
        const int64_t pixelJ(pixelLogicalIndex.getJ());
        m_image->setPixelColor(pixelI, pixelJ, QColor(pixelRGBA[0],
                                                      pixelRGBA[1],
                                                      pixelRGBA[2],
                                                      pixelRGBA[3]));
    }
}

/**
 * Set the pixel RGBA at the given pixel I and J.
 *
 * @param tabIndex
 *    Index of the tab.
 * @param overlayIndex
 *    Index of overlay
 * @param pixelIndex
 *     Pixel index with (0, 0) in top left corner
 * @param pixelRGBAOut
 *     RGBA at Pixel I, J
 */
void
ImageFile::setPixelRGBA(const int32_t /*tabIndex*/,
                        const int32_t /*overlayIndex*/,
                        const PixelIndex& pixelIndex,
                        const uint8_t pixelRGBA[4])
{
    if (m_image != NULL) {
        const int64_t pixelI(pixelIndex.getI());
        const int64_t pixelJ(pixelIndex.getJ());
        m_image->setPixelColor(pixelI, pixelJ, QColor(pixelRGBA[0],
                                                      pixelRGBA[1],
                                                      pixelRGBA[2],
                                                      pixelRGBA[3]));
    }
}

/**
 * Set the pixel RGBA at the given pixel I and J.
 *
 * @param tabIndex
 *    Index of the tab.
 * @param overlayIndex
 *    Index of overlay
 * @param pixelI
 *     Pixel I index (origin top left)
 * @param pixelJ
 *     Pixel J index (origin top left)
 * @param pixelRGBAOut
 *     RGBA at Pixel I, J
 */
void
ImageFile::setPixelRGBA(const int32_t /*tabIndex*/,
                        const int32_t /*overlayIndex*/,
                        const int32_t pixelI,
                        const int32_t pixelJ,
                        const uint8_t pixelRGBA[4])
{
    if (m_image != NULL) {
        m_image->setPixelColor(pixelI, pixelJ, QColor(pixelRGBA[0],
                                                      pixelRGBA[1],
                                                      pixelRGBA[2],
                                                      pixelRGBA[3]));
    }
}

/**
 * Set a full row of pixels with RGBA colors
 * @param rowIndex
 *    Index of the row
 * @param pixelRowRGBA
 *    Vector containing RGBA for all pixels in the row
 */
void
ImageFile::setPixelRowRGBA(const int32_t rowIndex,
                           const std::vector<uint8_t>& pixelRowRGBA)
{
    const int32_t rowLength(getWidth());
    if ((rowLength * 4) != static_cast<int32_t>(pixelRowRGBA.size())) {
        const AString msg("pixel RGBA vector is incorrect length="
                          + AString::number(pixelRowRGBA.size())
                          + ", should be="
                          + AString::number(rowLength * 4));
        CaretLogSevere(msg);
        CaretAssertMessage(0, msg);
        return;
    }
    if ((rowIndex < 0)
        || (rowIndex >= getHeight())) {
        const AString msg("Row index="
                          + AString::number(rowIndex)
                          + " out of range, height="
                          + AString::number(getHeight()));
        CaretLogSevere(msg);
        CaretAssertMessage(0, msg);
        return;
    }
    
    uint8_t* rowByteColor(m_image->scanLine(rowIndex));
    uint32_t* rowIntColor((uint32_t*)rowByteColor);
    switch (m_image->format()) {
        case QImage::Format_ARGB32:
            /*
             * This format is typically used by image in Workbench
             */
            for (int32_t i = 0; i < rowLength; i++) {
                const int32_t i4(i*4);
                rowIntColor[i] = QRgba64::fromRgba(pixelRowRGBA[i4],
                                                   pixelRowRGBA[i4+1],
                                                   pixelRowRGBA[i4+2],
                                                   pixelRowRGBA[i4+3]).toArgb32();
            }
            break;
        default:
        {
            const int32_t pixelJ(rowIndex);
            for (int32_t i = 0; i < rowLength; i++) {
                const int32_t i4(i*4);
                m_image->setPixelColor(i, pixelJ,
                                       QColor(pixelRowRGBA[i4],
                                              pixelRowRGBA[i4+1],
                                              pixelRowRGBA[i4+2],
                                              pixelRowRGBA[i4+3]));
            }
        }
            break;
    }
}

/**
 * Get the RGBA bytes from the image resized into the given width and height.
 *
 * @param imageOrigin
 *    Location of first pixel in the image data.
 * @param resizeToWidth
 *    New width for image.
 * @param resizeToHeight
 *    New height of the image.
 * @param bytesRGBAOut
 *    The RGBA bytes in the image.
 * @return
 *    True if the bytes, width, and height are valid, else false.
 */
bool
ImageFile::getImageResizedBytes(const IMAGE_DATA_ORIGIN_LOCATION imageOrigin,
                                const int32_t resizeToWidth,
                                const int32_t resizeToHeight,
                                std::vector<uint8_t>& bytesRGBAOut) const
{
    CaretAssertMessage(0, "There is a bug in this function, it has not been fixed. "
                       "Use getGraphicsPrimitiveForFeaturesImageDrawing()");
    
    bytesRGBAOut.clear();
    
    if (m_image == NULL) {
        return false;
    }
    
    const int32_t numBytes = resizeToWidth * resizeToHeight;
    if (numBytes <= 0) {
        return false;
    }
    
    const int32_t colorComponentsPerByte = 4;
    bytesRGBAOut.resize(numBytes * colorComponentsPerByte);
    
    QImage scaledImage = m_image->scaled(resizeToWidth,
                                   resizeToHeight,
                                   Qt::IgnoreAspectRatio,
                                   Qt::SmoothTransformation);
    CaretAssert(scaledImage.width() == resizeToWidth);
    CaretAssert(scaledImage.height() == resizeToHeight);
    
    /*
     * QImage::scaled() failed.
     */
    if (scaledImage.isNull()) {
        return false;
    }
    
    bool isOriginAtTop = false;
    switch (imageOrigin) {
        case IMAGE_DATA_ORIGIN_AT_BOTTOM:
            isOriginAtTop = false;
            break;
        case IMAGE_DATA_ORIGIN_AT_TOP:
            isOriginAtTop = true;
            break;
    }
    
    /*
     * Documentation for QImage states that setPixel may be very costly
     * and recommends using the scanLine() method to access pixel data.
     */
    for (int y = 0; y < resizeToHeight; y++) {
        const int scanLineIndex = (isOriginAtTop
                                   ? y
                                   : resizeToHeight -y - 1);
        const uchar* scanLine = scaledImage.scanLine(scanLineIndex);
        QRgb* rgbScanLine = (QRgb*)scanLine;
        
        const int32_t rowLength(resizeToWidth * 4);
        const int32_t rowOffset(rowLength * y);
        for (int x = 0; x < resizeToWidth; x++) {
            const int32_t contentOffset(rowOffset + (x * 4));
            QRgb& rgb = rgbScanLine[x];
            CaretAssertVectorIndex(bytesRGBAOut, contentOffset + 3);
            bytesRGBAOut[contentOffset] = static_cast<uint8_t>(qRed(rgb));
            bytesRGBAOut[contentOffset+1] = static_cast<uint8_t>(qGreen(rgb));
            bytesRGBAOut[contentOffset+2] = static_cast<uint8_t>(qBlue(rgb));
            bytesRGBAOut[contentOffset+3] = static_cast<uint8_t>(qAlpha(rgb));
        }
    }
    
    {
        QImage test(bytesRGBAOut.data(),
                    resizeToWidth,
                    resizeToHeight,
                    scaledImage.format());
        test.save("john.png");
    }
    return true;
}

/**
 * @return width of media file
 */
int32_t
ImageFile::getWidth() const
{
    int32_t w = 0;
    
    if (m_image != NULL) {
        w = m_image->width();
    }
    
    return w;
}

/**
 * @return height of media file
 */
int32_t
ImageFile::getHeight() const
{
    int32_t h = 0;
    
    if (m_image != NULL) {
        h = m_image->height();
    }
    
    return h;
}

/**
 * Essentially writes the image file to a byte array using the given format.
 *
 * @param byteArrayOut
 *    Byte array into which the image is written.
 * @param format
 *    Format for the image (jpg, ppm, etc.).
 */
void
ImageFile::getImageInByteArray(QByteArray& byteArrayOut,
                               const AString& format) const
{
    byteArrayOut.clear();
    
    if (m_image != NULL) {
        QBuffer buffer(&byteArrayOut);
        if ( ! buffer.open(QIODevice::WriteOnly)) {
            throw DataFileException(getFileName(),
                                    "PROGRAM ERROR: Unable to open byte array for output of image.");
        }
        
        bool successFlag = false;
        if (format.isEmpty()) {
            successFlag = m_image->save(&buffer);
        }
        else {
            successFlag = m_image->save(&buffer,
                                        format.toLatin1().data());
        }
        
        if ( ! successFlag) {
            throw DataFileException(getFileName(),
                                    "Failed to write image to byte array.  "
                                    + buffer.errorString());
        }
    }
}

/**
 * Essentially reads the image file from a byte array using the given format.
 *
 * @param byteArray
 *    Byte array from which the image is read.
 * @param format
 *    Format for the image (jpg, ppm, etc.) or empty if unknown.
 */
bool
ImageFile::setImageFromByteArray(const QByteArray& byteArray,
                                 const AString& format)
{
    bool successFlag = false;
    if (format.isEmpty()) {
        successFlag = m_image->loadFromData(byteArray);
    }
    else {
        successFlag = m_image->loadFromData(byteArray,
                                            format.toLatin1().data());
    }
    
    if ( ! successFlag) {
        CaretLogSevere(getFileName()
                       + " Failed to create image from byte array.");
    }

    return successFlag;
}

/**
 * Convert this image into a Volume File using the 
 * encapsulated Control Point File whose matrix 
 * must have been updated.
 *
 * @param colorMode
 *     Color mode for conversion.
 * @param errorMessageOut
 *     Contains error message if conversion fails.
 * @return
 *     Pointer to volume file or NULL if there is an error.
 *     Name of volume file is the name of the image but
 *     the file extension is changed to a volume file extension.
 *
 */
VolumeFile*
ImageFile::convertToVolumeFile(const CONVERT_TO_VOLUME_COLOR_MODE colorMode,
                               AString& errorMessageOut) const
{
    errorMessageOut.clear();
    
    std::vector<uint8_t> rgbaBytes;
    int32_t width = 0;
    int32_t height = 0;
    getImageBytesRGBA(ImageFile::IMAGE_DATA_ORIGIN_AT_BOTTOM,
                      rgbaBytes,
                      width,
                      height);
    if ((width <= 0)
        || (height <= 0)) {
        errorMessageOut = "Image width and/or height is invalid.";
        return NULL;
    }
    
    const Matrix4x4* transformationMatrix = m_controlPointFile->getLandmarkTransformationMatrix();
    
    float firstPixel[3] = { 0, 0, 0 };
    transformationMatrix->multiplyPoint3(firstPixel);
    if (imageDebugFlag) std::cout << "First pixel coord: " << AString::fromNumbers(firstPixel, 3, ",") << std::endl;
    
    float lastPixel[3] = { (float)(width - 1), (float)(height - 1), 0.0f };
    transformationMatrix->multiplyPoint3(lastPixel);
    if (imageDebugFlag) std::cout << "Last pixel coord: " << AString::fromNumbers(lastPixel, 3, ",") << std::endl;
    
    {
        float bl[3] = { 0.0, 0.0, 0.0 };
        transformationMatrix->multiplyPoint3(bl);
        ControlPoint3D bottomLeft(0, 0, 0, bl[0], bl[1], bl[2]);
        
        float br[3] = { (float)(width - 1.0), 0.0, 0.0 };
        transformationMatrix->multiplyPoint3(br);
        ControlPoint3D bottomRight(width - 1.0, 0.0, 0.0, br[0], br[1], br[2]);
        
        float tr[3] = { (float)(width - 1.0), (float)(height - 1.0), 0.0 };
        transformationMatrix->multiplyPoint3(tr);
        ControlPoint3D topRight((float)(width - 1.0), (float)(height - 1.0), 0.0, tr[0], tr[1], tr[2]);
        
        ControlPointFile volumeControlPointFile;
        volumeControlPointFile.addControlPoint(bottomLeft);
        volumeControlPointFile.addControlPoint(bottomRight);
        volumeControlPointFile.addControlPoint(topRight);
        
        if ( ! volumeControlPointFile.updateLandmarkTransformationMatrix(errorMessageOut)) {
            errorMessageOut.insert(0, "Volume Matrix: ");
            return NULL;
        }
    }
    
    std::vector<int64_t> dimensions;
    dimensions.push_back(width);  // I
    dimensions.push_back(height); // J
    dimensions.push_back(1); // K
    
    /*
     * Convert matrix4x4 to volume file vector of vectors.
     */
    std::vector<float> row1;
    std::vector<float> row2;
    std::vector<float> row3;
    std::vector<float> row4;
    
    for (int j = 0; j < 4; j++) {
        row1.push_back(transformationMatrix->getMatrixElement(0, j));
        row2.push_back(transformationMatrix->getMatrixElement(1, j));
        row3.push_back(transformationMatrix->getMatrixElement(2, j));
        row4.push_back(transformationMatrix->getMatrixElement(3, j));
    }
    std::vector<std::vector<float> > indexToSpace;
    indexToSpace.push_back(row1);
    indexToSpace.push_back(row2);
    indexToSpace.push_back(row3);
    indexToSpace.push_back(row4);
    
    int64_t numComponents = 1;
    SubvolumeAttributes::VolumeType whatType = SubvolumeAttributes::FUNCTIONAL;
    switch (colorMode) {
        case CONVERT_TO_VOLUME_COLOR_GRAYSCALE:
            break;
        case CONVERT_TO_VOLUME_COLOR_RGB:
            numComponents = 3;
            whatType = SubvolumeAttributes::RGB;
            break;
    }
    VolumeFile* volumeFile = new VolumeFile(dimensions,
                                            indexToSpace,
                                            numComponents,
                                            whatType);
    
    FileInformation fileInfo(getFileName());
    const AString volumeFileName = FileInformation::assembleFileComponents(fileInfo.getAbsolutePath(),
                                                                           fileInfo.getFileNameNoExtension(),
                                                                           DataFileTypeEnum::toFileExtension(DataFileTypeEnum::VOLUME));
    volumeFile->setFileName(volumeFileName);
    
    int32_t rgbaIndex = 0;
    const int64_t k = 0;
    const int64_t mapIndex = 0;
    for (int64_t j = 0; j < height; j++) {
        for (int64_t i = 0; i < width; i++) {
            switch (colorMode) {
                case CONVERT_TO_VOLUME_COLOR_GRAYSCALE:
                {
                    /*
                     * Luminosity conversion from GIMP
                     * http://docs.gimp.org/2.6/en/gimp-tool-desaturate.html
                     */
                    float intensity = ((rgbaBytes[rgbaIndex] * 0.21)
                                       + (rgbaBytes[rgbaIndex + 1] * 0.72)
                                       + (rgbaBytes[rgbaIndex + 2] * 0.07));
                    if (intensity > 255.0) intensity = 255.0;
                    else if (intensity < 0.0) intensity = 0.0;
                    
                    if (rgbaBytes[rgbaIndex + 3] <= 0.0) {
                        intensity = 0.0;
                    }
                    volumeFile->setValue(intensity, i, j, k, mapIndex, 0);
                    
                    rgbaIndex += 4;
                }
                    break;
                case CONVERT_TO_VOLUME_COLOR_RGB:
                {
                    CaretAssertVectorIndex(rgbaBytes, rgbaIndex);
                    volumeFile->setValue(rgbaBytes[rgbaIndex], i, j, k, mapIndex, 0);
                    CaretAssertVectorIndex(rgbaBytes, rgbaIndex);
                    volumeFile->setValue(rgbaBytes[rgbaIndex+1], i, j, k, mapIndex, 1);
                    CaretAssertVectorIndex(rgbaBytes, rgbaIndex);
                    volumeFile->setValue(rgbaBytes[rgbaIndex+2], i, j, k, mapIndex, 2);
                    if (numComponents == 4) {
                        CaretAssertVectorIndex(rgbaBytes, rgbaIndex);
                        volumeFile->setValue(rgbaBytes[rgbaIndex+3], i, j, k, mapIndex, 3);
                    }
                    rgbaIndex += 4;
                }
                    break;
            }
        }
    }
    
    switch (colorMode) {
        case CONVERT_TO_VOLUME_COLOR_GRAYSCALE:
        {
            PaletteColorMapping* pcm = volumeFile->getMapPaletteColorMapping(mapIndex);
            pcm->setSelectedPaletteToGrayInterpolated();
            pcm->setDisplayNegativeDataFlag(false);
            pcm->setDisplayZeroDataFlag(false);
            pcm->setDisplayPositiveDataFlag(true);
            pcm->setScaleMode(PaletteScaleModeEnum::MODE_AUTO_SCALE);
        }
            break;
        case CONVERT_TO_VOLUME_COLOR_RGB:
            break;
    }
    
    volumeFile->clearVoxelColoringForMap(mapIndex);
    volumeFile->updateScalarColoringForMap(mapIndex);
    
    return volumeFile;
}

/**
 * @return The control point file.
 */
ControlPointFile*
ImageFile::getControlPointFile()
{
    return m_controlPointFile;
}

/**
 * @return The control point file.
 */
const ControlPointFile*
ImageFile::getControlPointFile() const
{
    return m_controlPointFile;
}

GraphicsPrimitiveV3fT2f*
ImageFile::getGraphicsPrimitiveForFeaturesImageDrawing() const
{
    if (m_image == NULL) {
        return NULL;
    }
    
    if (m_featuresImageGraphicsPrimitive == NULL) {
        GraphicsPrimitiveV3fT2f* primitive(createGraphicsPrimitive());
        m_featuresImageGraphicsPrimitive.reset(primitive);
    }
    
    const int32_t pixelPrimitiveVertexStartIndex(0);
    const int32_t pixelPrimitiveVertexCount(4);
    m_featuresImageGraphicsPrimitive->setDrawArrayIndicesSubset(pixelPrimitiveVertexStartIndex,
                                                                pixelPrimitiveVertexCount);
    return m_featuresImageGraphicsPrimitive.get();
}

/**
 * @return The graphics primitive for drawing the image as a texture in media drawing model.
 * @param tabIndex
 *    Index of tab where image is drawn
 * @param overlayIndex
 *    Index of overlay
 */
GraphicsPrimitiveV3fT2f*
ImageFile::getGraphicsPrimitiveForMediaDrawing(const int32_t /*tabIndex*/,
                                               const int32_t /*overlayIndex*/) const
{
    if (m_image == NULL) {
        return NULL;
    }
    
    if (m_graphicsPrimitive == NULL) {
        GraphicsPrimitiveV3fT2f* primitive(createGraphicsPrimitive());
        m_graphicsPrimitive.reset(primitive);
    }
    CaretAssert(m_pixelPrimitiveVertexStartIndex >= 0);
    CaretAssert(m_pixelPrimitiveVertexCount > 0);
    m_graphicsPrimitive->setDrawArrayIndicesSubset(m_pixelPrimitiveVertexStartIndex,
                                                   m_pixelPrimitiveVertexCount);
    return m_graphicsPrimitive.get();
}

/**
 * @return A new graphics primitive for loaded data for both coordinate types
 */
GraphicsPrimitiveV3fT2f*
ImageFile::createGraphicsPrimitive() const
{
    /*
     * If image is too big for OpenGL texture limits, scale image to acceptable size
     */
    const int32_t maxTextureWidthHeight = GraphicsUtilitiesOpenGL::getTextureWidthHeightMaximumDimension();
    if (maxTextureWidthHeight > 0) {
        const int32_t excessWidth(m_image->width() - maxTextureWidthHeight);
        const int32_t excessHeight(m_image->height() - maxTextureWidthHeight);
        if ((excessWidth > 0)
            || (excessHeight > 0)) {
            if (excessWidth > excessHeight) {
                CaretLogWarning(getFileName()
                                + " is too big for texture.  Maximum width/height is: "
                                + AString::number(maxTextureWidthHeight)
                                + " Image Width: "
                                + AString::number(m_image->width())
                                + " Image Height: "
                                + AString::number(m_image->height()));
            }
        }
    }
    
    /*
     * Format must be RGB or ARGB for compatibility with OpenGL
     */
    verifyFormatCompatibleWithOpenGL();
    
    const std::array<float, 4> textureBorderColorRGBA { 0.0, 0.0, 0.0, 0.0 };
    
    GraphicsTextureSettings::PixelFormatType pixelFormat(GraphicsTextureSettings::PixelFormatType::BGRA);
    switch (m_image->format()) {
        case QImage::Format_RGB32:  /* Contains alpha that is always 255 */
            pixelFormat = GraphicsTextureSettings::PixelFormatType::BGRX;
            break;
        case QImage::Format_RGB888:
            pixelFormat = GraphicsTextureSettings::PixelFormatType::RGB;
            break;
        case QImage::Format_ARGB32:
            pixelFormat = GraphicsTextureSettings::PixelFormatType::BGRA;
            break;
        default:
            CaretAssertMessage(0, "Format not compatible with OpenGL");
            break;
    }
    
    /*
     * Filtering for matching image pixel to screen pixel
     */
    GraphicsTextureMagnificationFilterEnum::Enum magFilter(GraphicsTextureMagnificationFilterEnum::NEAREST);
    GraphicsTextureMinificationFilterEnum::Enum minFilter(GraphicsTextureMinificationFilterEnum::NEAREST);
    const bool smoothFlag(false);
    if (smoothFlag) {
        /*
         * This will break opacity on drawing volume slices over
         * histology slices that results in white edges appearing
         * when opacity is less than one.
         */
        magFilter = GraphicsTextureMagnificationFilterEnum::LINEAR;
        minFilter = GraphicsTextureMinificationFilterEnum::LINEAR_MIPMAP_LINEAR;
    }
    /*
     * Compress texture if image is large and compression is enabled
     */
    const GraphicsTextureSettings::CompressionType textureCompressionType(isImageTextureCompressed()
                                                                          ? GraphicsTextureSettings::CompressionType::ENABLED
                                                                          : GraphicsTextureSettings::CompressionType::DISABLED);
    GraphicsTextureSettings textureSettings(m_image->constBits(),
                                            m_image->width(),
                                            m_image->height(),
                                            1, /* slices */
                                            GraphicsTextureSettings::DimensionType::FLOAT_STR_2D,
                                            pixelFormat,
                                            GraphicsTextureSettings::PixelOrigin::TOP_LEFT,
                                            GraphicsTextureSettings::WrappingType::CLAMP,
                                            GraphicsTextureSettings::MipMappingType::ENABLED,
                                            textureCompressionType,
                                            magFilter,
                                            minFilter,
                                            textureBorderColorRGBA);
    GraphicsPrimitiveV3fT2f* primitive = GraphicsPrimitive::newPrimitiveV3fT2f(GraphicsPrimitive::PrimitiveType::OPENGL_TRIANGLE_STRIP,
                                                                               textureSettings);
    
    /*
     * Create a primitive for PIXEL coordinates
     *
     * Coordinates at EDGE of the pixels
     */
    const float minX = 0;
    const float maxX = getWidth();
    const float minY = 0;
    const float maxY = getHeight();
    
    /*
     * A Triangle Strip (consisting of two triangles) is used
     * for drawing the image.
     * The order of the vertices in the triangle strip is
     * Top Left, Bottom Left, Top Right, Bottom Right.
     * ORIGIN IS AT TOP LEFT
     */
    const float minTextureST(0.0);
    const float maxTextureST(1.0);
    m_pixelPrimitiveVertexStartIndex = primitive->getNumberOfVertices();
    primitive->addVertex(minX, minY, minTextureST, minTextureST);  /* Top Left */
    primitive->addVertex(minX, maxY, minTextureST, maxTextureST);  /* Bottom Left */
    primitive->addVertex(maxX, minY, maxTextureST, minTextureST);  /* Top Right */
    primitive->addVertex(maxX, maxY, maxTextureST, maxTextureST);  /* Bottom Right */
    m_pixelPrimitiveVertexCount = (primitive->getNumberOfVertices()
                                   - m_pixelPrimitiveVertexStartIndex);
    
    /*
     * Create a primitive for plane coordinates if available
     */
    if (isPlaneXyzSupported()) {
        /*
         * A Triangle Strip (consisting of two triangles) is used
         * for drawing the image.
         * The order of the vertices in the triangle strip is
         * Top Left, Bottom Left, Top Right, Bottom Right.
         * ORIGIN IS AT TOP LEFT
         */
        const float minTextureST(0.0);
        const float maxTextureST(1.0);
        const Vector3D coordinateTopLeft(getPlaneXyzTopLeft());
        const Vector3D coordinateTopRight(getPlaneXyzTopRight());
        const Vector3D coordinateBottomLeft(getPlaneXyzBottomLeft());
        const Vector3D coordinateBottomRight(getPlaneXyzBottomRight());
        m_planePrimitiveVertexStartIndex = primitive->getNumberOfVertices();
        primitive->addVertex(coordinateTopLeft[0],     coordinateTopLeft[1],     minTextureST, minTextureST);  /* Top Left */
        primitive->addVertex(coordinateBottomLeft[0],  coordinateBottomLeft[1],  minTextureST, maxTextureST);  /* Bottom Left */
        primitive->addVertex(coordinateTopRight[0],    coordinateTopRight[1],    maxTextureST, minTextureST);  /* Top Right */
        primitive->addVertex(coordinateBottomRight[0], coordinateBottomRight[1], maxTextureST, maxTextureST);  /* Bottom Right */
        m_planePrimitiveVertexCount = (primitive->getNumberOfVertices()
                                       - m_planePrimitiveVertexStartIndex);
    }
    
    return primitive;
}

/**
 * Verify and convert format so compatible with OpenGL
 */
void
ImageFile::verifyFormatCompatibleWithOpenGL() const
{
    /*
     * Using RGB888 (3 bytes per pixel) in place of RGB32 (4 bytes per pixel
     * with an alpha that is always 0xff) saves 25% in memory usage
     */
    if ((m_image->format() != QImage::Format_RGB888)
        && (m_image->format() != QImage::Format_ARGB32)) {
        if (m_image->hasAlphaChannel()) {
#if QT_VERSION >= QT_VERSION_CHECK(5, 13, 0)
            m_image->convertTo(QImage::Format_ARGB32);
#else
            *m_image = m_image->convertToFormat(QImage::Format_ARGB32);
#endif
        }
        else {
#if QT_VERSION >= QT_VERSION_CHECK(5, 13, 0)
            m_image->convertTo(QImage::Format_RGB888);
#else
            *m_image = m_image->convertToFormat(QImage::Format_RGB888);
#endif
        }
    }
}

/**
 * Save file data from the scene.  For subclasses that need to
 * save to a scene, this method should be overriden.  sceneClass
 * will be valid and any scene data should be added to it.
 *
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    restoring the scene.
 *
 * @param sceneClass
 *     sceneClass to which data members should be added.
 */
void
ImageFile::saveFileDataToScene(const SceneAttributes* sceneAttributes,
                                   SceneClass* sceneClass)
{
    MediaFile::saveFileDataToScene(sceneAttributes,
                                   sceneClass);
    
    if (m_controlPointFile != NULL) {
        sceneClass->addClass(m_controlPointFile->saveToScene(sceneAttributes,
                                                             "m_controlPointFile"));
    }
    
    /*
     * Added 01 June 2021 to assist with scenes created before default scaling
     */
    sceneClass->addInteger(ImageFile::SCENE_VERSION_NUMBER, 1);
}

/**
 * Restore file data from the scene.  For subclasses that need to
 * restore from a scene, this method should be overridden. The scene class
 * will be valid and any scene data may be obtained from it.
 *
 * @param sceneAttributes
 *    Attributes for the scene.  Scenes may be of different types
 *    (full, generic, etc) and the attributes should be checked when
 *    restoring the scene.
 *
 * @param sceneClass
 *     sceneClass for the instance of a class that implements
 *     this interface.  Will NEVER be NULL.
 */
void
ImageFile::restoreFileDataFromScene(const SceneAttributes* sceneAttributes,
                                        const SceneClass* sceneClass)
{
    MediaFile::restoreFileDataFromScene(sceneAttributes,
                                        sceneClass);
    
    m_controlPointFile->restoreFromScene(sceneAttributes,
                                         sceneClass->getClass("m_controlPointFile"));
    
    
    /*
     * If not restored, it will cause a warning
     *const int32_t sceneVersionNumber =
     */
    sceneClass->getIntegerValue(ImageFile::SCENE_VERSION_NUMBER, 0);
}

/**
 * @return File casted to an image file (avoids use of dynamic_cast that can be slow)
 */
ImageFile*
ImageFile::castToImageFile()
{
    return this;
}

/**
 * @return File casted to an image file (avoids use of dynamic_cast that can be slow)
 * Overidden in ImageFile
 */
const ImageFile*
ImageFile::castToImageFile() const
{
    return this;
}

/**
 * @param dataFileInformation
 *    Item to which information is added.
 */
void
ImageFile::addToDataFileContentInformation(DataFileContentInformation& dataFileInformation)
{
    MediaFile::addToDataFileContentInformation(dataFileInformation);
    
    if (m_image != NULL) {
        dataFileInformation.addNameAndValue("Width (pixels)", m_image->width());
        dataFileInformation.addNameAndValue("Height (pixels)", m_image->height());
        
        const float dotsPerMeter(m_image->dotsPerMeterX());
        if (dotsPerMeter > 0.0) {
            dataFileInformation.addNameAndValue("Width (meters)", m_image->width()   / dotsPerMeter);
            dataFileInformation.addNameAndValue("Height (meters)", m_image->height() / dotsPerMeter);
            
            /*
             * "To" and "From" units are flipped since conversion is on "per unit"
             */
            const float dotsPerInch = UnitsConversion::convertLength(UnitsConversion::LengthUnits::INCHES,
                                                                     UnitsConversion::LengthUnits::METERS,
                                                                     dotsPerMeter);
            dataFileInformation.addNameAndValue("Width (inches)", m_image->width()   / dotsPerInch);
            dataFileInformation.addNameAndValue("Height (inches)", m_image->height() / dotsPerInch);
            dataFileInformation.addNameAndValue("Pixels Per Meter", dotsPerMeter);
            dataFileInformation.addNameAndValue("Pixels Per Inch", dotsPerInch);
        }
        else {
            dataFileInformation.addNameAndValue("Pixels Per Meter", "Unavailable");
        }
        
        addPlaneCoordsToDataFileContentInformation(dataFileInformation);
        
        dataFileInformation.addNameAndValue("Color Table", (m_image->colorTable().empty()
                                                            ? "No"
                                                            : "Yes"));
    }
}

/**
 * Get all image file extensions supported by Qt for reading and writing image files.
 * @param readableExtensionsOut
 *    Output contains all readable image file extensions
 * @param writableExtensionsOut
 *    Output contains all writable image file extensions
 */
void
ImageFile::getQtSupportedImageFileExtensions(std::vector<AString>& readableExtensionsOut,
                                             std::vector<AString>& writableExtensionsOut)
{
    DataFileTypeEnum::getQtSupportedImageFileExtensions(readableExtensionsOut,
                                                        writableExtensionsOut);
}

/**
 * Get all image file extensions supported by Qt that support the 'clipRect' option
 * to read a portion of the image file.
 * @param clipRectReadableExtensionsOut
 *    Output contains all readable image file extensions that support 'clipRect'
 * @param scaledClipRectReadableExtensionsOut
 *    Output contains all readable image file extensions that support 'scaledClipRect'
 * @param metaDataReadableWritableExtensionsOut
 * Output contains all image file extensions that support  file metadata (the 'description' option) for reading and writing
 */
void
ImageFile::getImageFileQtSupportedOptionExtensions(std::vector<AString>& clipRectReadableExtensionsOut,
                                                   std::vector<AString>& scaledClipRectReadableExtensionsOut,
                                                   std::vector<AString>& metaDataReadableWritableExtensionsOut)
{
    clipRectReadableExtensionsOut.clear();
    scaledClipRectReadableExtensionsOut.clear();
    metaDataReadableWritableExtensionsOut.clear();
    
    std::vector<AString> readableExtensions;
    std::vector<AString> writableExtensions;
    DataFileTypeEnum::getQtSupportedImageFileExtensions(readableExtensions,
                                                        writableExtensions);
    
    for (auto& ext : writableExtensions) {
        if (std::find(readableExtensions.begin(),
                      readableExtensions.end(),
                      ext) != readableExtensions.end()) {
            QString name(QDir::tempPath() + "/file." + ext);
            
            /*
             * Image file must exist for testing supported options
             */
            QImage image(2, 2, QImage::Format_RGB32);
            QImageWriter writer(name);
            if (writer.canWrite()) {
                writer.write(image);
                
                QImageReader reader(name);
                if (reader.supportsOption(QImageIOHandler::ClipRect)) {
                    clipRectReadableExtensionsOut.push_back(ext);
                }
                if (reader.supportsOption(QImageIOHandler::ScaledClipRect)) {
                    scaledClipRectReadableExtensionsOut.push_back(ext);
                }
                if (reader.supportsOption(QImageIOHandler::Description)
                    && writer.supportsOption(QImageIOHandler::Description)) {
                    metaDataReadableWritableExtensionsOut.push_back(ext);
                }
            }
        }
    }
}

/**
 * Get all image file extensions supported by Workbench for reading and writing image files.
 * These are a subset of the extensions supported by Qt.
 * @param readableExtensionsOut
 *    Output contains all readable image file extensions
 * @param writableExtensionsOut
 *    Output contains all writable image file extensions
 * @param defaultWritableExtension
 *    The default file extension
 */
void
ImageFile::getWorkbenchSupportedImageFileExtensions(std::vector<AString>& readableExtensionsOut,
                                                    std::vector<AString>& writableExtensionsOut,
                                                    AString& defaultWritableExtension)
{
    DataFileTypeEnum::getWorkbenchSupportedImageFileExtensions(readableExtensionsOut,
                                                               writableExtensionsOut,
                                                               defaultWritableExtension);
}

/**
 * Move the QImage metadata to the file metadata
 */
void
ImageFile::readFileMetaDataFromQImage()
{
    CaretAssert(m_fileMetaData);
    m_fileMetaData->clear();
    
    if (m_image != NULL) {
        if ( ! m_image->isNull()) {
            Matrix4x4 scaledToPlaneMatrix;
            bool scaledToPlaneMatrixValidFlag(false);
            Matrix4x4 planeToMillimetersMatrix;
            bool planeToMillimetersMatrixValidFlag(false);
            
            QStringList allKeys = m_image->textKeys();
            QStringListIterator iter(allKeys);
            while (iter.hasNext()) {
                const AString key(iter.next());
                const AString value(m_image->text(key));
                m_fileMetaData->set(key, value);
                
                if (key == getMetaDataNameScaledToPlaneMatrix()) {
                    scaledToPlaneMatrix.setMatrixFromRowMajorOrderString(value);
                    scaledToPlaneMatrixValidFlag = true;
                }
                else if (key == getMetaDataNamePlaneToMillimetersMatrix()) {
                    planeToMillimetersMatrix.setMatrixFromRowMajorOrderString(value);
                    planeToMillimetersMatrixValidFlag = true;
                }
            }
            
            if (scaledToPlaneMatrixValidFlag
                && planeToMillimetersMatrixValidFlag) {
                std::shared_ptr<CziNonLinearTransform> unusedOne, unusedTwo;
                setTransformMatrices(scaledToPlaneMatrix,
                                     scaledToPlaneMatrixValidFlag,
                                     planeToMillimetersMatrix,
                                     planeToMillimetersMatrixValidFlag,
                                     unusedOne,
                                     unusedTwo);
            }
        }
    }
    
    m_fileMetaData->clearModified();
}

/**
 * Move the file metadata into the QImage
 */
void
ImageFile::writeFileMetaDataToQImage() const
{
    if (m_image != NULL) {
        /*
         * May clear metadata, but not sure
         */
        m_image->setText("", "");
        
        const auto allKeys(m_fileMetaData->getAllMetaDataNames());
        for (const auto& key : allKeys) {
            const AString value = m_fileMetaData->get(key);
            m_image->setText(key, value);
        }
    }
}

/**
 * @return True if this file supports file metadata, else false.
 * Support of metadata is dependent upon the file's extension.  Only
 * some image file types support metadata.
 */
bool
ImageFile::supportsFileMetaData() const
{
    std::vector<AString> clipRectReadableExtensions;
    std::vector<AString> scaledClipRectReadableExtensions;
    std::vector<AString> metaDataReadableWritableExtensions;
    ImageFile::getImageFileQtSupportedOptionExtensions(clipRectReadableExtensions,
                                                       scaledClipRectReadableExtensions,
                                                       metaDataReadableWritableExtensions);

    FileInformation fileInfo(getFileName());
    const AString ext(fileInfo.getFileExtension());
    
    if (std::find(metaDataReadableWritableExtensions.begin(),
                  metaDataReadableWritableExtensions.end(),
                  ext) != metaDataReadableWritableExtensions.end()) {
        return true;
    }
    
    return false;
}

/**
 * @return True if the given pixel index is valid for the CZI image file (may be outside of currently loaded sub-image)
 * @param frameIndex
 *    Index of frame
 * @param pixelIndexOriginAtTopLeft
 *    Image of pixel with origin (0, 0) at the top left
 */
bool
ImageFile::isPixelIndexInFrameValid(const int32_t /*frameIndex*/,
                                    const PixelIndex& pixelIndexOriginAtTopLeft) const
{
    return isPixelIndexValid(pixelIndexOriginAtTopLeft);
}

/**
 * @return True if the given pixel index is valid
 * @param frameIndex
 *    Index of frame
 * @param pixelLogicalIndex
 *    Pixel logical index
 */
bool
ImageFile::isPixelIndexInFrameValid(const int32_t /*frameIndex*/,
                                    const PixelLogicalIndex& pixelLogicalIndex) const
{
    return isPixelIndexValid(pixelLogicalIndex);
}


/**
 * @return True if the given pixel index is valid for the image in the given tabf
 * @param pixelIndexOriginAtTopLeft
 *    Image of pixel with origin (0, 0) at the top left
 */
bool
ImageFile::isPixelIndexValid(const PixelIndex& pixelIndexOriginAtTopLeft) const
{
    if (m_image != NULL) {
        const int32_t i(pixelIndexOriginAtTopLeft.getI());
        const int32_t j(pixelIndexOriginAtTopLeft.getJ());
        if ((i >= 0)
            && (i < m_image->width())
            && (j >= 0)
            && (j < m_image->height())) {
            return true;
        }
    }
    return false;
}

/**
 * @return True if the given pixel index is valid
 * @param pixelLogicalIndex
 *    Pixel logical index
 */
bool
ImageFile::isPixelIndexValid(const PixelLogicalIndex& pixelLogicalIndex) const
{
    if (m_image != NULL) {
        const int32_t i(pixelLogicalIndex.getI());
        const int32_t j(pixelLogicalIndex.getJ());
        if ((i >= 0)
            && (i < m_image->width())
            && (j >= 0)
            && (j < m_image->height())) {
            return true;
        }
    }
    return false;
}

/**
 * Get the identification text for the pixel at the given pixel index with origin at bottom left.
 * @param tabIndex
 *    Index of the tab in which identification took place
 * @param frameIndices
 *    Indics of the frames
 * @param pixelLogicalIndex
 *    Logical pixel index
 * @param columnOneTextOut
 *    Text for column one that is displayed to user.
 * @param columnTwoTextOut
 *    Text for column two that is displayed to user.
 * @param toolTipTextOut
 *    Text for tooltip
 */
void
ImageFile::getPixelLogicalIdentificationTextForFrames(const int32_t tabIndex,
                                               const std::vector<int32_t>& frameIndices,
                                                  const PixelLogicalIndex& pixelLogicalIndex,
                                                  std::vector<AString>& columnOneTextOut,
                                                  std::vector<AString>& columnTwoTextOut,
                                                  std::vector<AString>& toolTipTextOut) const
{
    columnOneTextOut.clear();
    columnTwoTextOut.clear();
    toolTipTextOut.clear();
    
    std::vector<int32_t> validFrameIndices;
    for (int32_t frameIndex : frameIndices) {
        if (isPixelIndexInFrameValid(frameIndex,
                                     pixelLogicalIndex)) {
            validFrameIndices.push_back(frameIndex);
        }
    }
    if (validFrameIndices.empty()) {
        return;
    }
    
    columnOneTextOut.push_back("Filename");
    columnTwoTextOut.push_back(getFileNameNoPath());
    
    const AString pixelText("Pixel IJ ("
                            + AString::number(pixelLogicalIndex.getI())
                            + ","
                            + AString::number(pixelLogicalIndex.getJ())
                            + ")");
    toolTipTextOut.push_back(pixelText);
    
    const float dotsPerMeter(m_image->dotsPerMeterX());
    if (dotsPerMeter > 0.0) {
        /*
         * Show coordinates in millimeters
         */
        const float metersPerDot(1 / dotsPerMeter);
        const float xMillimeters(pixelLogicalIndex.getI() * metersPerDot * 100.0);
        const float yMillimeters(pixelLogicalIndex.getJ() * metersPerDot * 100.0);
        const AString mmText("("
                             + AString::number(xMillimeters, 'f', 2)
                             + "mm, "
                             + AString::number(yMillimeters, 'f', 2)
                             + "mm)");
        columnOneTextOut.push_back(pixelText);
        columnTwoTextOut.push_back(mmText);
        
        toolTipTextOut.push_back(mmText);
    }
    else {
        columnOneTextOut.push_back(pixelText);
        columnTwoTextOut.push_back("");
    }
    
    for (int32_t frameIndex : validFrameIndices) {
        uint8_t rgba[4];
        if (getPixelRGBA(tabIndex,
                         frameIndex,
                         pixelLogicalIndex,
                         rgba)) {
            const AString leftText("Frame "
                                   + AString::number(frameIndex + 1));
            const AString rightText("RGBA ("
                                    + AString::fromNumbers(rgba, 4, ",")
                                    + ")");
            columnOneTextOut.push_back(leftText);
            columnTwoTextOut.push_back(rightText);
            toolTipTextOut.push_back(leftText
                                     + " "
                                     + rightText);
        }
    }

    CaretAssert(columnOneTextOut.size() == columnTwoTextOut.size());
}

/**
 * Find the Pixel nearest the given XYZ coordinate
 * @param xyz
 *    The coordinate
 * @param signedDistanceToPixelMillimetersOut
 *    Output with signed distance to the pixel in millimeters
 * @param pixelLogicalIndexOut
 *    Output with logical pixel index
 * @return
 *    True if successful, else false.
 */
bool
ImageFile::findPixelNearestStereotaxicXYZ(const Vector3D& xyz,
                                          float& signedDistanceToPixelMillimetersOut,
                                          PixelLogicalIndex& pixelLogicalIndexOut) const
{
    return MediaFile::findPixelNearestStereotaxicXYZ(xyz,
                                                     signedDistanceToPixelMillimetersOut,
                                                     pixelLogicalIndexOut);
}

/*
 * @return Primitive for drawing media with coordinates
 * @param tabIndex
 *    Index of tab where image is drawn
 * @param overlayIndex
 *    Index of overlay
 */
GraphicsPrimitiveV3fT2f*
ImageFile::getGraphicsPrimitiveForPlaneXyzDrawing(const int32_t /*tabIndex*/,
                                                  const int32_t /*overlayIndex*/) const
{
    if (m_image == NULL) {
        return NULL;
    }
    
    if ( ! isPlaneXyzSupported()) {
        return NULL;
    }
    
    if (m_graphicsPrimitive == NULL) {
        GraphicsPrimitiveV3fT2f* primitive(createGraphicsPrimitive());
        m_graphicsPrimitive.reset(primitive);
    }
    CaretAssert(m_planePrimitiveVertexStartIndex >= 0);
    CaretAssert(m_planePrimitiveVertexCount > 0);
    m_graphicsPrimitive->setDrawArrayIndicesSubset(m_planePrimitiveVertexStartIndex,
                                                   m_planePrimitiveVertexCount);
    return m_graphicsPrimitive.get();
}

/**
 * @return True if the texture for the image should be compressed to save memory.
 */
bool
ImageFile::isImageTextureCompressed() const
{
    EventCaretPreferencesGet prefsEvent;
    EventManager::get()->sendEvent(prefsEvent.getPointer());
    CaretPreferences* prefs = prefsEvent.getCaretPreferences();
    CaretAssert(prefs);
    
    if (prefs->isImageFileTextureCompressionEnabled()) {
        const int64_t megabyte(1000000);
        const int64_t byteSize(getTextureCompressionSizeMegabytes() * megabyte);
        const int64_t imageBytes(getWidth() * getHeight() * 4);
        if (imageBytes > byteSize) {
            return true;
        }
    }
    return false;
}

