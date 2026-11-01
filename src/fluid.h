// SPDX-FileCopyrightText: 2026 ficool2
// SPDX-License-Identifier: MIT

#pragma once

typedef struct b3World b3World;

// Apply buoyancy and drag from all fluid sensors
void b3ApplyFluidForces( b3World* world, float timeStep );

// Wake the sleeping bodies in fluids that changed or are moving
void b3WakeFluidVisitors( b3World* world );
