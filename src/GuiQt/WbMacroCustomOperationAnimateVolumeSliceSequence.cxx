
/*LICENSE_START*/
/*
 *  Copyright (C) 2019 Washington University School of Medicine
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

#define __WB_MACRO_CUSTOM_OPERATION_ANIMATE_VOLUME_SLICE_DECLARE__
#include "WbMacroCustomOperationAnimateVolumeSliceSequence.h"
#undef __WB_MACRO_CUSTOM_OPERATION_ANIMATE_VOLUME_SLICE_DECLARE__

#include "BrainBrowserWindow.h"
#include "BrowserTabContent.h"
#include "CaretAssert.h"
#include "Matrix4x4.h"
#include "Model.h"
#include "ModelVolume.h"
#include "ModelWholeBrain.h"
#include "Overlay.h"
#include "OverlaySet.h"
#include "SystemUtilities.h"
#include "VolumeMappableInterface.h"
#include "WbMacroCustomDataTypeEnum.h"
#include "WbMacroCustomOperationTypeEnum.h"
#include "WuQMacroCommand.h"
#include "WuQMacroCommandParameter.h"
#include "WuQMacroExecutorMonitor.h"

using namespace caret;


    
/**
 * \class caret::WbMacroCustomOperationAnimateVolumeSliceSequence
 * \brief Macro custom operation incrementing volume slices
 * \ingroup GuiQt
 */

/**
 * Constructor.
 */
WbMacroCustomOperationAnimateVolumeSliceSequence::WbMacroCustomOperationAnimateVolumeSliceSequence()
: WbMacroCustomOperationBase(WbMacroCustomOperationTypeEnum::ANIMATE_VOLUME_SLICE_SEQUENCE)
{
    
}

/**
 * Destructor.
 */
WbMacroCustomOperationAnimateVolumeSliceSequence::~WbMacroCustomOperationAnimateVolumeSliceSequence()
{
}

/**
 * Get a new instance of the macro command
 *
 * @return
 *     Pointer to command or NULL if not valid
 *     Use getErrorMessage() for error information if NULL returned
 */
WuQMacroCommand*
WbMacroCustomOperationAnimateVolumeSliceSequence::createCommand()
{
    const int32_t versionOne(1);
    
    
    const QString description("Sequence through the volume slices: \n"
                              "(1) start at the \"selected slice\"; \n"
                              "(2) decrement the slice index to the first slice; \n"
                              "(3) increment the slice index to the last slice; \n"
                              "(4) decrement the slice index returning to the \"selected slice\"");
    QString errorMessage;
    WuQMacroCommand* command = WuQMacroCommand::newInstanceCustomCommand(WbMacroCustomOperationTypeEnum::toName(getOperationType()),
                                                                         versionOne,
                                                                         "none",
                                                                         WbMacroCustomOperationTypeEnum::toGuiName(getOperationType()),
                                                                         description,
                                                                         1.0,
                                                                         errorMessage);
    if (command != NULL) {
        WuQMacroCommandParameter* paramOne = new WuQMacroCommandParameter(WuQMacroDataValueTypeEnum::AXIS,
                                                                          "Volume Axis",
                                                                          "Z");
        command->addParameter(paramOne);
        command->addParameter(WuQMacroDataValueTypeEnum::FLOAT,
                              "Duration (secs)",
                              (float)20.0);
    }
    else {
        appendToErrorMessage(errorMessage);
    }
    
    return command;
}

/**
 * Execute the macro command
 *
 * @param parent
 *     Parent widget for any dialogs
 * @param executorMonitor
 *     the macro executor monitor
 * @param executorOptions
 *     Options for executor
 * @param macroCommand
 *     macro command to run
 * @return
 *     True if command executed successfully, else false
 *     Use getErrorMessage() for error information if false returned
 */
bool
WbMacroCustomOperationAnimateVolumeSliceSequence::executeCommand(QWidget* parent,
                                                           const WuQMacroExecutorMonitor* executorMonitor,
                                                           const WuQMacroExecutorOptions* executorOptions,
                                                           const WuQMacroCommand* macroCommand)
{
    CaretAssert(parent);
    CaretAssert(macroCommand);
    
    if ( ! validateCorrectNumberOfParameters(macroCommand, 2)) {
        return false;
    }
    const QString axisName(macroCommand->getParameterAtIndex(0)->getValue().toString().toUpper());
    const float durationSeconds = macroCommand->getParameterAtIndex(1)->getValue().toFloat();

    BrainBrowserWindow* bbw = qobject_cast<BrainBrowserWindow*>(parent);
    if (bbw == NULL) {
        appendToErrorMessage("Parent for running surface macro is not a browser window.");
        return false;
    }
    
    BrowserTabContent* tabContent = bbw->getBrowserTabContent();
    if (tabContent == NULL) {
        appendToErrorMessage("No tab is selected in browser window.");
        return false;
    }


    Axis axis = Axis::X;
    if (axisName == "X") {
        axis = Axis::X;
    }
    else if (axisName == "Y") {
        axis = Axis::Y;
    }
    else if (axisName == "Z") {
        axis = Axis::Z;
    }
    else {
        appendToErrorMessage("Axis named \""
                             + axisName
                             + "\" is invalid.  Use X, Y, or Z.");
    }
    
    if (durationSeconds < 0.0) {
        appendToErrorMessage("Duration must be greater than zero.");
    }
    
    Model* model = tabContent->getModelForDisplay();
    if (model == NULL) {
        appendToErrorMessage("No model for surface rotation");
    }

    if ( ! getErrorMessage().isEmpty()) {
        return false;
    }

    const bool successFlag = performSliceIncrement(executorMonitor,
                                                   executorOptions,
                                                   tabContent,
                                                   axis,
                                                   durationSeconds);
    return successFlag;
}

/**
 * @param executorMonitor
 *     The macro executor's monitor
 * @param executorOptions
 *     The executor options
 * @param tabContent
 *     Content in the selected tab
 * @param axis
 *     Axis for viewing
 * @param durationSecondes
 *     Duration of time for command to run
 */
bool
WbMacroCustomOperationAnimateVolumeSliceSequence::performSliceIncrement(const WuQMacroExecutorMonitor* executorMonitor,
                                                                  const WuQMacroExecutorOptions* executorOptions,
                                                                  BrowserTabContent* tabContent,
                                                                  const Axis axis,
                                                                  const float durationSeconds)
{
    ModelVolume* volumeModel(tabContent->getDisplayedVolumeModel());
    ModelWholeBrain* wholeBrainModel(tabContent->getDisplayedWholeBrainModel());
    if ((volumeModel == NULL)
        && (wholeBrainModel == NULL)) {
        appendToErrorMessage("For slice increment, View must be All or Volume");
        return false;
    }
    
    const VolumeMappableInterface* vmi = tabContent->getOverlaySet()->getUnderlayVolume();
    if (vmi == NULL) {
        appendToErrorMessage("No volume is selected as an overlay");
        return false;
    }
    std::vector<int64_t> dims;
    vmi->getDimensions(dims);
    if (dims.size() < 3) {
        appendToErrorMessage("Dimensions are invalid for underlay volume");
    }
    
    int32_t numberOfSlices(0);
    int32_t startingSliceIndex(0);
    switch (axis) {
        case Axis::X:
            numberOfSlices = dims[0];
            startingSliceIndex = tabContent->getVolumeSliceIndexParasagittal(vmi);
            break;
        case Axis::Y:
            numberOfSlices = dims[1];
            startingSliceIndex = tabContent->getVolumeSliceIndexCoronal(vmi);
            break;
        case Axis::Z:
            numberOfSlices = dims[2];
            startingSliceIndex = tabContent->getVolumeSliceIndexAxial(vmi);
            break;
    }
    if (numberOfSlices <= 0) {
        appendToErrorMessage("No slices in underlay volume for selected dimension");
    }
    
    const float defaultNumberOfSteps(numberOfSlices);
    float numberOfSteps(0.0);
    float iterationSleepTime(0.0);
    getNumberOfStepsAndSleepTime(executorOptions,
                                 defaultNumberOfSteps,
                                 durationSeconds,
                                 numberOfSteps,
                                 iterationSleepTime);

    enum SliceMode {
        DECREMENT_TO_ZERO,
        INCREMENT_TO_LAST,
        DECREMENT_TO_START
    };
    SliceMode sliceMode = DECREMENT_TO_ZERO;
    
    float sliceIndex = startingSliceIndex;
    float sliceIncrement = numberOfSlices / (numberOfSteps / 2);
    
    bool doneFlag(false);
    while ( ! doneFlag) {
        switch (axis) {
            case Axis::X:
                tabContent->setVolumeSliceIndexParasagittal(vmi, sliceIndex);
                break;
            case Axis::Y:
                tabContent->setVolumeSliceIndexCoronal(vmi, sliceIndex);
                break;
            case Axis::Z:
                tabContent->setVolumeSliceIndexAxial(vmi, sliceIndex);
                break;
        }

        switch (sliceMode) {
            case DECREMENT_TO_START:
                sliceIndex -= sliceIncrement;
                if (sliceIndex <= startingSliceIndex) {
                    doneFlag = true;
                }
                break;
            case DECREMENT_TO_ZERO:
                sliceIndex -= sliceIncrement;
                if (sliceIndex <= 0.0) {
                    sliceIndex = 0;
                    sliceMode = INCREMENT_TO_LAST;
                }
                break;
            case INCREMENT_TO_LAST:
                sliceIndex += sliceIncrement;
                if (sliceIndex >= (numberOfSlices - 1)) {
                    sliceIndex = (numberOfSlices - 1);
                    sliceMode = DECREMENT_TO_START;
                }
                break;
        }
        

        updateGraphics();        
        
        if (executorMonitor->testForStop()) {
            appendToErrorMessage(executorMonitor->getStoppedByUserMessage());
            return false;
        }
        
        sleepForSecondsAtEndOfIteration(iterationSleepTime);
    }
    
    switch (axis) {
        case Axis::X:
            tabContent->setVolumeSliceIndexParasagittal(vmi, startingSliceIndex);
            break;
        case Axis::Y:
            tabContent->setVolumeSliceIndexCoronal(vmi, startingSliceIndex);
            break;
        case Axis::Z:
            tabContent->setVolumeSliceIndexAxial(vmi, startingSliceIndex);
            break;
    }
    updateGraphics();
    
    return true;
}
