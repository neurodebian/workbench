
/*LICENSE_START*/
/*
 *  Copyright (C) 2021 Washington University School of Medicine
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

#define __GRAPHICS_FRAMES_PER_SECOND_DECLARE__
#include "GraphicsFramesPerSecond.h"
#undef __GRAPHICS_FRAMES_PER_SECOND_DECLARE__

#include <algorithm>

#include "CaretAssert.h"
#include "ElapsedTimer.h"

using namespace caret;


    
/**
 * \class caret::GraphicsFramesPerSecond 
 * \brief Assists with timing graphics
 * \ingroup Graphics
 *
 * There are three frame-per-second timers in an instance of this class
 *
 * (1) A timer that reports the last drawing time (time between calls to startOfDrawing()
 * and endOfDrawing().
 * (2) A timer that reports an average drawing time (average of time between calls to
 * startOfDrawing() and endOfDrawing().  Call getAverageFramesPerSecond().
 * (3) A timer that reports the average time between calls to endOfDrawing().
 * Call getAverageIntervalFramesPerSecond().
 *
 */

/**
 * Constructor.
 * @param maximumFrameCount
 *    Number of frames that are averaged to create the frames per second
 */
GraphicsFramesPerSecond::GraphicsFramesPerSecond(const int32_t maximumFrameCount)
: CaretObject(),
m_maximumFrameCount(maximumFrameCount)
{
    m_startEndTimer.reset(new ElapsedTimer());
    m_intervalTimer.reset(new ElapsedTimer());
    m_intervalFrameTimes.resize(m_maximumFrameCount);
    std::fill(m_intervalFrameTimes.begin(),
              m_intervalFrameTimes.end(),
              -1.0);
    
    reinitialize();
}

/**
 * Destructor.
 */
GraphicsFramesPerSecond::~GraphicsFramesPerSecond()
{
}

/**
 * @return The most recent frames per second between calls to 'startOfDrawing()' and
 * 'endOfDrawing()'.
 */
double
GraphicsFramesPerSecond::getFramesPerSecond() const
{
    return m_mostRecentFramesPerSecond;
}

/**
 * @return The average frames per second between calls to 'startOfDrawing()' and
 * 'endOfDrawing()'.
 */
double
GraphicsFramesPerSecond::getAverageFramesPerSecond() const
{
    return getFPS(m_startStopFrameTimes);
}

/**
 * @return Average of time since between calls to 'endOfDrawing()'.
 */
double
GraphicsFramesPerSecond::getAverageIntervalFramesPerSecond() const
{
    double fpsOut(0.0);
    
    if (m_intervalTimer->isStarted()) {
        fpsOut = getFPS(m_intervalFrameTimes);
        
    }
    
    return fpsOut;
}

/**
 * Compute the frames per second for a given group of frame times.
 * @param frameTimes
 *    The frame times
 * @return The frames per second
 */
double
GraphicsFramesPerSecond::getFPS(const std::vector<double>& frameTimes) const
{
    /*
     * Get sum of recent frame times
     */
    double sumMilliseconds(0.0);
    int32_t numValidFrames(0);
    for (int32_t i = 0; i < m_maximumFrameCount; i++) {
        CaretAssertVectorIndex(frameTimes, i);
        const double& ft = frameTimes[i];
        if (ft >= 0.0) {
            sumMilliseconds += ft;
            numValidFrames++;
        }
    }
    
    /*
     * Calculate frames per second
     */
    double framesPerSecond(0.0);
    const double seconds(sumMilliseconds / 1000.0);
    if (seconds > 0.0) {
        framesPerSecond = (numValidFrames / seconds);
    }
    
    return framesPerSecond;
}

/**
 * Call at the beginning of drawing to start timing one frame
 */
void
GraphicsFramesPerSecond::startOfDrawing()
{
    /*
     * Start timing
     */
    m_startEndTimer->reset();
    m_startEndTimer->start();
}

/**
 * Call at the end of drawing to end timing for one frame
 */
void
GraphicsFramesPerSecond::endOfDrawing()
{
    const double elapsedTime(m_startEndTimer->getElapsedTimeMilliseconds());
    
    /*
     * Save frame time
     */
    if (m_startStopFrameTimesIndex >= m_maximumFrameCount) {
        m_startStopFrameTimesIndex = 0;
    }
    CaretAssertVectorIndex(m_startStopFrameTimes, m_startStopFrameTimesIndex);
    m_startStopFrameTimes[m_startStopFrameTimesIndex] = elapsedTime;
    ++m_startStopFrameTimesIndex;
    
    /*
     * Convert most recent to frames per second
     */
    m_mostRecentFramesPerSecond = 0.0;
    if (elapsedTime > 0.0) {
        const double seconds = (elapsedTime / 1000.0);
        m_mostRecentFramesPerSecond = 1.0 / seconds;
    }
    
    updateAverageIntervalFramesPerSecond();
}

/**
 * Reinitialize the instance
 */
void
GraphicsFramesPerSecond::reinitialize()
{
    m_startStopFrameTimes.resize(m_maximumFrameCount);
    std::fill(m_startStopFrameTimes.begin(), m_startStopFrameTimes.end(), -1.0);
    m_startStopFrameTimesIndex = 0;
    m_mostRecentFramesPerSecond = 0.0;
}

/**
 * Get a description of this object's content.
 * @return String describing this object's content.
 */
AString 
GraphicsFramesPerSecond::toString() const
{
    return "GraphicsFramesPerSecond";
}

/**
 * Update average of time since this function was last called.
 */
void
GraphicsFramesPerSecond::updateAverageIntervalFramesPerSecond()
{
    if (m_intervalTimer->isStarted()) {
        /*
         * Since graphics are updated (as needed), we exclude any
         * long updates that are most likley due to:
         * (1) The user interacting (rotate, pan, etc);
         * (2) User does nothing for a while;
         * (3) User starts interacting again.
         *
         * If we do not exclude these long updates we will
         * have some big jumps in the timing since the time
         * doing nothing would be measured
         */
        const double userDoesNothingMilliseconds(5.0 * 1000.0);
        const double intervalMilliseconds(m_intervalTimer->getElapsedTimeMilliseconds());
        if (intervalMilliseconds < userDoesNothingMilliseconds) {
            if (m_intervalFrameTimesIndex >= m_maximumFrameCount) {
                m_intervalFrameTimesIndex = 0;
            }
            CaretAssertVectorIndex(m_intervalFrameTimes, m_intervalFrameTimesIndex);
            m_intervalFrameTimes[m_intervalFrameTimesIndex] = intervalMilliseconds;
            ++m_intervalFrameTimesIndex;
        }
    }

    m_intervalTimer->reset();
    m_intervalTimer->start();
}

