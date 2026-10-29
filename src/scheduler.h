// SPDX-FileCopyrightText: 2026 Erin Catto
// SPDX-License-Identifier: MIT

#pragma once

#include "box3d/box3d.h"

// Starts a step, i.e. resets the task slots.
void b3BeginSchedulerStep( b3Scheduler* scheduler );
void b3EndSchedulerStep( b3Scheduler* scheduler );

// See b3EnqueueTaskCallback and b3FinishTaskCallback
void* b3SchedulerEnqueueTask( b3TaskCallback* task, void* taskContext, void* userContext, const char* name );
void b3SchedulerFinishTask( void* userTask, void* userContext );
