
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
#include <iostream>

#define __CIFTI_CONNECTIVITY_MATRIX_DENSE_DYNAMIC_FILE_DECLARE__
#include "CiftiConnectivityMatrixDenseDynamicFile.h"
#undef __CIFTI_CONNECTIVITY_MATRIX_DENSE_DYNAMIC_FILE_DECLARE__

#include "CaretAssert.h"
#include "CaretLogger.h"
#include "CaretOMP.h"
#include "CiftiBrainordinateDataSeriesFile.h"
#include "CiftiFile.h"
#include "ElapsedTimer.h"
#include "FileInformation.h"

using namespace caret;


    
/**
 * \class caret::CiftiConnectivityMatrixDenseDynamicFile 
 * \brief Connectivity Dynamic Dense x Dense File version of data-series
 * \ingroup Files
 *
 * Contains dynamic connectivity from brainordinates to brainordinates.
 * Internally, the file format is the same as a data series file.  When
 * a row is requested, the row is correlated with all other rows
 * producing the connectivity from that row to all other rows.
 */

/**
 * Constructor.
 *
 * @param parentDataSeriesFile
 *     Parent data series file.
 */
CiftiConnectivityMatrixDenseDynamicFile::CiftiConnectivityMatrixDenseDynamicFile(const CiftiBrainordinateDataSeriesFile* parentDataSeriesFile)
: CiftiMappableConnectivityMatrixDataFile(DataFileTypeEnum::CONNECTIVITY_DENSE_DYNAMIC),
m_parentDataSeriesFile(parentDataSeriesFile),
m_parentDataSeriesCiftiFile(NULL),
m_numberOfBrainordinates(-1),
m_numberOfTimePoints(-1),
m_enabledForUser(true),
m_cacheDataFlag(false)
{
    CaretAssert(m_parentDataSeriesFile);
}

/**
 * Destructor.
 */
CiftiConnectivityMatrixDenseDynamicFile::~CiftiConnectivityMatrixDenseDynamicFile()
{
    
}

/**
 * @return True if this file type supports writing, else false.
 *
 * Dense files do NOT support writing.
 */
bool
CiftiConnectivityMatrixDenseDynamicFile::supportsWriting() const
{
    return false;
}

/**
 * Is the file enabled for the user so that file is in user-interface.
 */
bool
CiftiConnectivityMatrixDenseDynamicFile::isEnabledForUser() const
{
    return m_enabledForUser;
}

/**
 * Set the file enabled for the user so that file is in user-interface.
 *
 * @param enabled
 *     File enabled status.
 */
void
CiftiConnectivityMatrixDenseDynamicFile::setEnabledForUser(const bool enabled)
{
    m_enabledForUser = enabled;
}


/**
 * This method is intended for overriding by subclasess so that they
 * can examine and verify the data that was read.  This method is
 * called after successfully reading a file.
 */
void
CiftiConnectivityMatrixDenseDynamicFile::updateAfterReading(const CiftiFile* ciftiFile)
{
    m_parentDataSeriesCiftiFile = const_cast<CiftiFile*>(ciftiFile);
    
    AString path, nameNoExt, ext;
    FileInformation fileInfo(m_parentDataSeriesCiftiFile->getFileName());
    fileInfo.getFileComponents(path, nameNoExt, ext);
    setFileName(FileInformation::assembleFileComponents(path,
                                                        nameNoExt,
                                                        DataFileTypeEnum::toFileExtension(DataFileTypeEnum::CONNECTIVITY_DENSE_DYNAMIC)));
    
    /*
     * Need dimensions of data
     * Note that CIFTI XML in this file is identifical to CIFTI XML in parent data-series file
     */
    const CiftiXML& ciftiXML = getCiftiFile()->getCiftiXML();
    m_numberOfBrainordinates = ciftiXML.getBrainModelsMap(CiftiXML::ALONG_COLUMN).getLength();
    m_numberOfTimePoints     = ciftiXML.getSeriesMap(CiftiXML::ALONG_ROW).getLength();
    
    m_rowData.clear();
    
    if ((m_numberOfBrainordinates > 0)
        && (m_numberOfTimePoints > 0)) {
        m_rowData.resize(m_numberOfBrainordinates);
        
        if (m_cacheDataFlag) {
            /*
             * Read all of the data.  Time-series type files are not
             * too large and by caching the data, it eliminates
             * numerous calls to read the data when correlation
             * is performed.
             *
             * READ DATA FROM PARENT FILE
             */
            for (int32_t i = 0; i < m_numberOfBrainordinates; i++) {
                CaretAssertVectorIndex(m_rowData, i);
                m_rowData[i].m_data.resize(m_numberOfTimePoints);
                m_parentDataSeriesCiftiFile->getRow(&m_rowData[i].m_data[0],
                                                    i);
            }
        }
        
        preComputeRowMeanAndSumSquared();
    }
}


/**
 * Load data for the given column.
 *
 * @param dataOut
 *     Output with data.
 * @param index of the column.
 */
void
CiftiConnectivityMatrixDenseDynamicFile::getDataForColumn(float* /*dataOut*/, const int64_t& /*index*/) const
{
    const AString msg("Should never be called for Dense Dynamic File");
    CaretAssertMessage(0, msg);
    CaretLogSevere(msg);
}

/**
 * Load data for the given row.
 *
 * @param dataOut
 *     Output with data.
 * @param index of the row.
 */
void
CiftiConnectivityMatrixDenseDynamicFile::getDataForRow(float* dataOut, const int64_t& index) const
{
    if ((m_numberOfBrainordinates <= 0)
        || (m_numberOfTimePoints <= 0)) {
        return;
    }
    
    ElapsedTimer timer;
    timer.start();
    
// disable for timing     #pragma omp CARET_PARFOR
    for (int32_t iRow = 0; iRow < m_numberOfBrainordinates; iRow++) {
        float coefficient = 1.0;
        
        if (iRow != index) {
            coefficient = correlation(index, iRow, m_numberOfTimePoints);
        }
        
        dataOut[iRow] = coefficient;
    }
    
    std::cout << "Time to correlate (milliseconds): " << timer.getElapsedTimeMilliseconds() << std::endl;
}

/**
 * Compute the mean and sum-squared for each row so that they
 * are only calculated once.
 */
void
CiftiConnectivityMatrixDenseDynamicFile::preComputeRowMeanAndSumSquared()
{
    CaretAssert(m_numberOfBrainordinates > 0);
    CaretAssert(m_numberOfTimePoints > 0);
    
    const float numPointsFloat = m_numberOfTimePoints;
    
    for (int32_t iRow = 0; iRow < m_numberOfBrainordinates; iRow++) {
        double sum = 0.0;
        double sumSquared = 0.0;

        CaretAssertVectorIndex(m_rowData, iRow);
        if (m_cacheDataFlag) {
            for (int32_t iPoint = 0; iPoint < m_numberOfTimePoints; iPoint++) {
                CaretAssertVectorIndex(m_rowData[iRow].m_data, iPoint);
                const float d = m_rowData[iRow].m_data[iPoint];
                sum        += d;
                sumSquared += (d * d);
            }
        }
        else {
            std::vector<float> data(m_numberOfTimePoints);
            m_parentDataSeriesCiftiFile->getRow(&data[0], iRow);
            for (int32_t iPoint = 0; iPoint < m_numberOfTimePoints; iPoint++) {
                CaretAssertVectorIndex(data, iPoint);
                const float d = data[iPoint];
                sum        += d;
                sumSquared += (d * d);
            }
        }
        
        
        const float mean       = (sum / numPointsFloat);
        m_rowData[iRow].m_mean = mean;
        const float ssxx = (sumSquared - (numPointsFloat * mean * mean));
        CaretAssert(ssxx >= 0.0);
        m_rowData[iRow].m_sqrt_ssxx = std::sqrt(ssxx);
    }
}

/**
 * Correlation from https://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
 *
 * @param rowIndex
 *     Index of a row
 * @param y
 *     Index of another row
 * @param numberOfPoints
 *     Number of points int the two arrays
 * @return
 *     The correlation coefficient computed on the two arrays.
 */
float
CiftiConnectivityMatrixDenseDynamicFile::correlation(const int32_t rowIndex,
                                                     const int32_t otherRowIndex,
                                                     const int32_t numberOfPoints) const
{
    const double numFloat = numberOfPoints;
    double xySum = 0.0;
    
    CaretAssertVectorIndex(m_rowData, rowIndex);
    CaretAssertVectorIndex(m_rowData, otherRowIndex);
    const RowData& data = m_rowData[rowIndex];
    const RowData& otherData = m_rowData[otherRowIndex];
    
    if (m_cacheDataFlag) {
        for (int i = 0; i < numberOfPoints; i++) {
            CaretAssertVectorIndex(data.m_data, i);
            CaretAssertVectorIndex(otherData.m_data, i);
            xySum += data.m_data[i] * otherData.m_data[i];
        }
    }
    else {
        std::vector<float> dataVector(m_numberOfTimePoints);
        std::vector<float> otherDataVector(m_numberOfTimePoints);
        m_parentDataSeriesCiftiFile->getRow(&dataVector[0], rowIndex);
        m_parentDataSeriesCiftiFile->getRow(&otherDataVector[0], otherRowIndex);
        
        for (int i = 0; i < numberOfPoints; i++) {
            CaretAssertVectorIndex(dataVector, i);
            CaretAssertVectorIndex(otherDataVector, i);
            xySum += dataVector[i] * otherDataVector[i];
        }
    }
    
    const double ssxy = xySum - (numFloat * data.m_mean * otherData.m_mean);
    
    float correlationCoefficient = 0.0;
    if ((data.m_sqrt_ssxx > 0.0)
        && (otherData.m_sqrt_ssxx > 0.0)) {
        correlationCoefficient = (ssxy / (data.m_sqrt_ssxx * otherData.m_sqrt_ssxx));
    }
    return correlationCoefficient;
}

