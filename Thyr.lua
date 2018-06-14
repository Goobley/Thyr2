--- A tool for simulating 3D gyrosynchrotron emission.
-- C M J Osborne - University of Glasgow, 2018
--
-- MIT License
--
-- @module Thyr


-- Libraries
local List = require('pl.List')
local pretty = require('pl.pretty')
local copy = require('pl.tablex').copy
local deepcopy = require('pl.tablex').deepcopy
local Plyght = require('Plyght.Plyght')
local maf = require('maf')
local threads = require('threads')

-- require 'jit'
-- jit.off()

-- Import C prototypes and ffi library
local ffi = require('ffi')
local defs = assert(io.open('Simulation.i', 'r'))
local defsStr = defs:read('*all')
defs:close()
ffi.cdef(defsStr)
ffi.cdef[[
void* calloc(size_t num, size_t size);
void free(void* p);
]]
local GyroIn = ffi.typeof('GyroSimDataC')
local ThermalRadiationData = ffi.typeof('ThermalRadiationData')
local DoublePtr = ffi.typeof('f64*') 
local Double = ffi.typeof('f64') 
local Gyro = ffi.load('./libgyro.so')

-- Useful constants.
local ArcToCm = math.pi / 180.0 / 3600.0 * 1.49597870e13
local SolarRadiusCm = 960 * ArcToCm
local MultiSampleFactor = 2
math.inf = math.huge

--- Allocates the 3D Voxel Grid to hold j and k.
-- Wrapper around the C functions that allocate the grid holding the
-- emission/absorption parameters, so that grids are automatically garbage
-- collected when there are no remaining references.
-- @int xSize Number of x voxels
-- @int ySize Number of y voxels
-- @int zSize Number of z voxels
-- @return Garbage collected grid
local function allocate_grid(xSize, ySize, zSize)
    return ffi.gc(Gyro.AllocateGrid(xSize, ySize, zSize), Gyro.FreeGrid)
end

--- GC'd array of Doubles.
-- Lua "type" for a simple garbage collected array of doubles. _Much_ smaller
-- than using a table.
-- @int len length of array (0-indexed as it comes from C)
local function DoubleCArray(len)
    return ffi.gc(ffi.cast(DoublePtr, ffi.C.calloc(len, ffi.sizeof(Double))), ffi.C.free)
end

--- Generator to index a uniform 3D grid from x, y, and z.
-- The grids are one-dimensional for speed reasons, this hides that fact.
-- Also does some error chacking to avoid segfaults.
-- @int numVox grid side length
-- @return function (x, y, z) -> index for use with uniform 3D grid
local function uniform_3d_index_gen(numVox)
    local f = function(x, y, z)
        assert(x < numVox, 'Check translation vector')
        assert(y < numVox, 'Check translation vector')
        assert(z < numVox, 'Check translation vector')
        return (z) * numVox^2 + (y) * numVox + x
    end
    return f 
end

--- Generator to index a 3D grid from x, y, and z.
-- The grids are one-dimensional for speed reasons, this hides that fact.
-- In this version x, y and z can have different lengths.
-- Also does some error chacking to avoid segfaults.
-- @int xVox grid x length
-- @int yVox grid y length
-- @int zVox grid z length
-- @return function (x, y, z) -> index for use with 3D grid
local function non_uniform_3d_index_gen(xVox, yVox, zVox)
    local f = function(x, y, z)
        assert(x < xVox, 'Check translation vector')
        assert(y < yVox, 'Check translation vector')
        assert(z < zVox, 'Check translation vector')
        if x < 0 then print('x', x) end
        if y < 0 then print('y', y) end
        if z < 0 then print('z', z) end
        local val =  (z) * xVox * yVox + (y) * xVox + x
        if val < 0 then print(x, y , z) assert(false, 'borked') end
        return val
    end
    return f 
end

--- Generator to index a 2D grid from x and y.
-- The grids are one-dimensional for speed reasons, this hides that fact.
-- This is used for things like images, (assumed to be square in Thyr).
-- @int len grid side length
-- @return function (x, y) -> index for use with uniform 2D grid
local function uniform_2d_index_gen(len)
    local f = function(x, y)
        return (y-1) * len + x
    end
    return f
end

--- Computes size to allocate from an AABB object and a scale factor.
-- @tparam AABB aabb an axis-aligned bounding box object (see later)
-- @int scale the scale factor by which the linear grid density is scaled
-- @return number of voxels to allocate
local function aabb_to_grid_alloc_size(aabb, scale)
    local xRange = math.ceil(scale * (aabb.x.max - aabb.x.min))
    local yRange = math.ceil(scale * (aabb.y.max - aabb.y.min))
    local zRange = math.ceil(scale * (aabb.z.max - aabb.z.min))

    return (xRange+1) * (yRange+1) * (zRange+1)
end

--- Index function from an AABB object and a scale factor.
-- @tparam AABB aabb an axis-aligned bounding box object (see later)
-- @int scale the scale factor by which the linear grid density is scaled
-- @return Index function for this grid
local function aabb_to_index(aabb, scale)
    local xRange = math.ceil(scale * (aabb.x.max - aabb.x.min))
    local yRange = math.ceil(scale * (aabb.y.max - aabb.y.min))
    local zRange = math.ceil(scale * (aabb.z.max - aabb.z.min))

    return non_uniform_3d_index_gen(xRange, yRange, zRange)
end

--- Serialise grid to file.
-- This saves two files determined by the filename prefix
-- These can be restored using @{restore_backup_grid}
-- @tparam Grid grid the grid object generated by @{dipole_in_aabb}
-- @tparam GyroDataTable gyroParams the shared data required for Gyro calculations
-- @tparam string filename a filename prefix, '_grid.bak' and '_meta.bak' will be appended for the files that are saved
local function backup_grid(grid, gyroParams, filename)
    local result = Gyro.SerializeGrid(filename..'_grid.bak', grid.grid, aabb_to_grid_alloc_size(grid.aabb, grid.scale), #gyroParams.frequency)
    if result ~= 0 then
        assert(false, 'Unable to serialize grid')
    end

    -- local gp = pretty.write(gyroParams);
    local g = deepcopy(grid)
    g.grid = nil
    g.idx = nil
    g.trans = {g.trans.x, g.trans.y, g.trans.z}
    g.rayDir = {g.rayDir.x, g.rayDir.y, g.rayDir.z}
    g.rotMat = {g.rotMat.x, g.rotMat.y, g.rotMat.z, g.rotMat.w}

    local ser = pretty.write({ grid = g, gyroParams = gyroParams })

    local f = assert(io.open(filename..'_meta.bak', 'w'))
    f:write(ser)
    f:flush()
    f:close()
end

--- Load grid from file.
-- Loads the two files determined by the filename prefix
-- @tparam string filename a filename prefix, '_grid.bak' and '_meta.bak' will be appended to load the files as per @{backup_grid}
-- @treturn Grid loaded grid
-- @treturn GyroParamTable the shared data required for Gyro calculations
local function restore_backup_grid(filename)
    local f = assert(io.open(filename..'_meta.bak', 'r'))
    local metaStr = f:read('*all')
    f:close()

    local meta = pretty.read(metaStr)
    local gyroParams = meta.gyroParams
    local grid = meta.grid
    grid.trans = maf.vec3(unpack(grid.trans))
    grid.rayDir = maf.vec3(unpack(grid.rayDir))
    grid.rotMat = maf.quat(unpack(grid.rotMat))
    -- VoxToCm = grid.params.gridSide / grid.params.numVox


    grid.grid = Gyro.DeserializeGrid(filename..'_grid.bak', aabb_to_grid_alloc_size(grid.aabb, grid.scale), #gyroParams.frequency)
    assert(grid.grid ~= nil, 'Something went wrong deserializing the grid')
    grid.idx = aabb_to_index(grid.aabb, grid.scale)

    return grid, gyroParams
end

--- Generate a dipole shape with the given parameters.
-- This shape is stored in an array of bytes with a 1 signaling that the voxel
-- is inside the dipole shape
-- @tparam DipoleParams d parameters determining shape of Dipole. See Driver function for form
-- @tparam[opt] maf.vec3 trans translation of dipole inside grid
-- @treturn ByteArray Grid of data describing dipole
-- @treturn IndexFn Function to use for indexing the ByteArray
-- @treturn AABB AABB for dipole contained in array
-- @treturn int Number of voxels that are inside the dipole
local function dipole(d, trans)
    trans = trans and trans or maf.vec3()

    local grid = ffi.new("uint8_t[?]", (d.numVox+1)^3)

    local sin, cos = math.sin, math.cos
    local min, max = math.min, math.max

    local index = uniform_3d_index_gen(d.numVox)

    local xRange = {min = math.inf, max = -math.inf}
    local yRange = {min = math.inf, max = -math.inf}
    local zRange = {min = math.inf, max = -math.inf}
    local xiRange = {min = math.inf, max = -math.inf}
    local etaRange = {min = math.inf, max = -math.inf}
    local SolarRadiusVox = SolarRadiusCm / d.VoxToCm

    local step = 0.5
    local count = 0
    for xi = - d.numVox/2, d.numVox/2, step do
        for eta = 0, d.numVox, step do
            for zeta = -d.numVox/2, d.numVox/2, step do
                
                local dipoleCond = (d.height - eta * (1.0 + (xi^2 / (zeta^2 + eta^2)))^1.5)^2
                                        + zeta^2 * (1.0 + (xi^2 / (zeta^2 + eta^2)))^3

                local x = xi * cos(d.phi0) - eta * sin(d.phi0)
                local y = xi * sin(d.phi0) + eta * cos(d.phi0)
                local z = zeta

                local dipoleCond2 = x^2 + (y - d.depth + SolarRadiusVox)^2 + z^2

                if dipoleCond <= d.rho0^2 and dipoleCond2 > SolarRadiusVox^2 then
                    local xx = math.floor(x + d.numVox / 2 + trans.x)
                    local yy = math.floor(y + trans.y)
                    local zz = math.floor(z + d.numVox / 2 + trans.z)
                    grid[index(xx, yy, zz)] = 1
                    count = count + 1
                    xRange.min = min(xRange.min, xx)
                    xRange.max = max(xRange.max, xx)
                    yRange.min = min(yRange.min, yy)
                    yRange.max = max(yRange.max, yy)
                    zRange.min = min(zRange.min, zz)
                    zRange.max = max(zRange.max, zz)
                    xiRange.min = min(xiRange.min, xi)
                    xiRange.max = max(xiRange.max, xi)
                    etaRange.min = min(etaRange.min, eta)
                    etaRange.max = max(etaRange.max, eta)
                end
            end
        end
    end
    local eps = 2
    xRange.min = xRange.min - eps
    xRange.max = xRange.max + eps
    yRange.min = yRange.min - eps
    yRange.max = yRange.max + eps
    zRange.min = zRange.min - eps
    zRange.max = zRange.max + eps

    return grid, index, {x = xRange, y = yRange, z = zRange, xi = xiRange, eta = etaRange}, count
end

--- Fill the dipole shape with the value of the requested parameter.
-- @tparam Grid d table containing the dipole shape and the parameters
-- @int scale the scale factor by which to increase the linear grid density
-- @tparam GyroParamTable gyroParams shared parameters used in Gyro calculations
-- @tparam string paramName parameter name to fill the grid with. One of 'angle', 'height', 'N_el', 'N_p', 'B_mag', 'temperature', 'N_HI', 'N_HII' 
-- @tparam function(height)->PlasmaParams interp_fn function returning, for a given height, the values of N_el, N_p, temperature, N_HI, N_HII
-- @treturn Grid a complete grid structure with the filled parameter 
local function visualise_dipole_param(d, scale, gyroParams, paramName, interp_fn)
    local p = d.params
    d.trans = d.trans and d.trans or maf.vec3()

    local xRange = math.ceil(scale * (d.aabb.x.max - d.aabb.x.min))
    local yRange = math.ceil(scale * (d.aabb.y.max - d.aabb.y.min))
    local zRange = math.ceil(scale * (d.aabb.z.max - d.aabb.z.min))

    local grid = DoubleCArray((xRange+1) * (yRange+1) * (zRange+1))

    local rSun = 960
    local sin, cos = math.sin, math.cos
    local min, max = math.min, math.max

    local gyroIn = GyroIn()
    gyroIn.aniso.type = 'None' 
    local thermalData = ThermalRadiationData()

    local index = non_uniform_3d_index_gen(xRange, yRange, zRange)
    local SolarRadiusVox = SolarRadiusCm / p.VoxToCm

    local step = 0.5 / scale
    local arghDex = 1
    -- for xiIdx = aabb.x.min, aabb.x.max, step do
    io.write('\n')
    for x = d.aabb.x.min - d.trans.x, d.aabb.x.max - d.trans.x, step do
        -- Use ANSI escape code. Don't be fooled by the \033 you see in C examples, C
        -- takes these escaped characters in octal, lua is in decimal 33octal = 27dec
        io.write('\027[A\027[K')
        io.write("Approximate: ", arghDex, "/", p.count * scale^3, "\n")
        -- for eta = aabb.y.min, aabb.y.max, step do
        for y = d.aabb.y.min - d.trans.y, d.aabb.y.max - d.trans.y, step do
            for zetaIdx = d.aabb.z.min - d.trans.z, d.aabb.z.max - d.trans.z, step do
                
                local xi = (x - p.numVox / 2) * cos(p.phi0) + y * sin(p.phi0)
                local eta = y * cos(p.phi0) - (x - p.numVox / 2) * sin(p.phi0)
                local zeta = zetaIdx - p.numVox / 2

                local dipoleCond = (p.height - eta * (1.0 + (xi^2 / (zeta^2 + eta^2)))^1.5)^2
                                        + zeta^2 * (1.0 + (xi^2 / (zeta^2 + eta^2)))^3

                -- local x = xi * cos(phi0) - eta * sin(phi0)
                -- local y = xi * sin(phi0) + eta * cos(phi0)
                local z = zetaIdx

                local dipoleCond2 = (x - p.numVox/ 2)^2 + (y - p.depth + SolarRadiusVox)^2 + (z - p.numVox/ 2)^2

                if dipoleCond <= p.rho0^2 and dipoleCond2 > SolarRadiusVox^2 then
                    -- -- Altitude of a cell is its y coordinate * vox2cm
                    -- -- The issue with this is that the AABB is lower than the start of the grid...
                    -- -- So these values start too high...
                    -- -- local h = (y - d.aabb.y.min + d.trans.y) * p.VoxToCm
                    -- -- Let's start by calculating the y coordinate of the surface for this x and z
                    --     local X = x - p.numVox / 2
                    --     local Z = z - p.numVox / 2
                    --     local b = 2 * (SolarRadiusVox - p.depth)
                    --     local Delta = b^2 - 4 * (p.depth^2 - 2 * p.depth * SolarRadiusVox + X^2 + Z^2)
                    -- local surfaceY = (-b + math.sqrt(Delta)) / 2

                    -- We already have a depth parameter... Duh
                    local h = (y - p.depth + d.trans.y) * p.VoxToCm
                    -- local h = (y - surfaceY) * p.VoxToCm
                    local xx = math.floor(scale * x - scale * d.aabb.x.min + scale * d.trans.x)
                    local yy = math.floor(scale * y - scale * d.aabb.y.min + scale * d.trans.y)
                    local zz = math.floor(scale * z - scale * d.aabb.z.min + scale * d.trans.z)

                    local b0x = cos(p.phi0) * p.b0Mag
                    local b0y = sin(p.phi0) * p.b0Mag

                    local xAxis = x - p.numVox / 2
                    local yAxis = y
                    local zAxis = z - p.numVox / 2

                    local b0 = maf.vec3(b0x, b0y, 0)
                    local rVec = maf.vec3(xAxis, yAxis, zAxis)
                    local mu = -b0 * p.height^3
                    local b = (rVec * (3 * mu:dot(rVec)) - mu * rVec:length()^2) / (rVec:length()^5)
                    local bMag = b:length()
                    gyroIn.bMag = bMag

                    local bDir = b:normalize():dot(d.rayDir)
                    local bAngle = math.acos(bDir)
                    gyroIn.angle = bAngle * 180 / math.pi

                    if gyroIn.angle > 87 and gyroIn.angle < 90 then
                        gyroIn.angle = 87
                    elseif gyroIn.angle >= 90 and gyroIn.angle < 93 then
                        gyroIn.angle = 93
                    elseif gyroIn.angle > 177 and gyroIn.angle < 180 then
                        gyroIn.angle = 177
                    elseif gyroIn.angle >= 180 and gyroIn.angle < 183 then
                        gyroIn.angle = 183
                    elseif gyroIn.angle < 3 then
                        gyroIn.angle = 3
                    end

                    local data = interp_fn(h)
                    gyroIn.nel = data.nel and data.nel or 0
                    gyroIn.np = data.np
                    thermalData.temperature = data.temperature
                    thermalData.protonDensity = data.HII
                    thermalData.neutralHDensity = data.HI

                    if (x > d.aabb.x.min and x < d.aabb.x.max) and
                       (y > d.aabb.y.min and y < d.aabb.y.max) and
                       (z > d.aabb.y.min and z < d.aabb.z.max)
                    then
                        local param = {
                                       ['angle'] = gyroIn.angle,
                                       ['height'] = h,
                                       ['N_el'] = gyroIn.nel,
                                       ['N_p'] = gyroIn.np,
                                       ['B_mag'] = gyroIn.bMag,
                                       ['N_HI'] = thermalData.neutralHDensity,
                                       ['temperature'] = thermalData.temperature,
                                       ['N_HII'] = thermalData.protonDensity
                                      }
                        grid[index(xx, yy, zz)] = param[paramName]
                    end
                    arghDex = arghDex + 1
                end
            end
        end
    end

    return {
            grid = grid,
            aabb = d.aabb,
            idx = index,
            scale = scale,
            params = p,
            rayDir = d.rayDir,
            trans = d.trans,
            rotMat = d.rotMat
        } 
end

--- Construct a thread pool to parallelise computation of the gyrosynchrotron coefficients
-- @int number of threads to use
-- @treturn threads.Threads Thread pool object used in torch with initialised threads
local function create_gyro_thread_pool(numThreads)
    print('Spawning '..numThreads..' threads for gyrosynchrotron calculations')
    return threads.Threads(numThreads,
                    function()
                        local ffi = require('ffi')
                        C = {}
                        ffi.cdef(defsStr)
                        ffi.cdef[[
                        void* calloc(size_t num, size_t size);
                        void free(void* p);
                        ]]
                        C.GyroIn = ffi.typeof('GyroSimDataC')
                        C.ThermalRadiationData = ffi.typeof('ThermalRadiationData')
                        C.DoublePtr = ffi.typeof('f64*') 
                        C.Double = ffi.typeof('f64') 
                        C.Gyro = ffi.load('./libgyro.so')
                        C.ffi = ffi
                    end)
end

--- Construct a dipole with j and k in the AABB.
-- @tparam Grid d the Grid returned from @{dipole}
-- @int scale the scale factor to apply
-- @tparam GyroParamTable gyroParams the shared parameters for the Gyro calculation
-- @tparam function(height)->PlasmaParams interp_fn function returning, for a given height, the values of N_el, N_p, temperature, N_HI, N_HII
-- @treturn Grid a complete grid structure with the filled j and k values for each requested frequency for each voxel within the dipole, these are null pointers elsewhere
local function dipole_in_aabb(d, scale, gyroParams, interp_fn, pool)
    -- collectgarbage('stop')
    local p = d.params
    d.trans = d.trans and d.trans or maf.vec3()

    local xRange = math.ceil(scale * (d.aabb.x.max - d.aabb.x.min))
    local yRange = math.ceil(scale * (d.aabb.y.max - d.aabb.y.min))
    local zRange = math.ceil(scale * (d.aabb.z.max - d.aabb.z.min))

    local grid = allocate_grid(xRange+1, yRange+1, zRange+1)
    -- local sampleCount = torch.IntTensor(xRange, yRange, zRange)
    local sampleCount = torch.IntTensor((xRange+1) * (yRange+1) * (zRange+1)):zero()
    local emptySamples = torch.IntTensor((xRange+1) * (yRange+1) * (zRange+1)):zero()

    local rSun = 960
    local sin, cos = math.sin, math.cos
    local min, max = math.min, math.max

    -- local gyroIn = GyroIn()
    -- local delta = DoubleCArray(#gyroParams.delta)
    -- gyroIn.delta = delta
    -- for i = 0,#gyroParams.delta-1 do delta[i] = gyroParams.delta[i+1] end
    -- gyroIn.deltaLen = #gyroParams.delta

    -- local delta = torch.Tensor(gyroParams.delta)
    -- local energy = torch.Tensor(gyroParams.energy)
    local frequency = torch.Tensor(gyroParams.frequency)


    -- local energy = DoubleCArray(#gyroParams.energy)
    -- gyroIn.energy = energy
    -- for i = 0,#gyroParams.energy-1 do energy[i] = gyroParams.energy[i+1] end
    -- gyroIn.energyLen = #gyroParams.energy

    -- local frequency = DoubleCArray(#gyroParams.frequency)
    -- gyroIn.frequency = frequency
    -- for i = 0,#gyroParams.frequency-1 do frequency[i] = gyroParams.frequency[i+1] end
    -- gyroIn.frequencyLen = #gyroParams.frequency

    -- local thermalData = ThermalRadiationData()

    local index = non_uniform_3d_index_gen(xRange, yRange, zRange)
    local SolarRadiusVox = SolarRadiusCm / p.VoxToCm

    local step = 1.0 / MultiSampleFactor / scale
    local arghDex = 1
    io.write('\n')
    local freeList = List {}

    for x = d.aabb.x.min - d.trans.x + 0.0*step, d.aabb.x.max - d.trans.x, step do
        io.write('\027[A\027[K')
        io.write("Approximate: ", arghDex, "/", p.count * scale^3, "\n")
        -- for eta = aabb.y.min, aabb.y.max, step do
        for y = d.aabb.y.min - d.trans.y + 0.0*step, d.aabb.y.max - d.trans.y, step do
            for zetaIdx = d.aabb.z.min - d.trans.z + 0.0*step, d.aabb.z.max - d.trans.z, step do
                
                local xi = (x - p.numVox / 2) * cos(p.phi0) + y * sin(p.phi0)
                local eta = y * cos(p.phi0) - (x - p.numVox / 2) * sin(p.phi0)
                local zeta = zetaIdx - p.numVox / 2

                local dipoleCond = (p.height - eta * (1.0 + (xi^2 / (zeta^2 + eta^2)))^1.5)^2
                                        + zeta^2 * (1.0 + (xi^2 / (zeta^2 + eta^2)))^3

                -- local x = xi * cos(phi0) - eta * sin(phi0)
                -- local y = xi * sin(phi0) + eta * cos(phi0)
                local z = zetaIdx

                local dipoleCond2 = (x - p.numVox/ 2)^2 + (y - p.depth + SolarRadiusVox)^2 + (z - p.numVox/ 2)^2

                if dipoleCond <= p.rho0^2 and dipoleCond2 > SolarRadiusVox^2 then
                    -- Altitude of a cell is its y coordinate * vox2cm
                    -- local h = (y - d.aabb.y.min + d.trans.y) * p.VoxToCm
                    local h = (y - p.depth + d.trans.y) * p.VoxToCm
                    local data = interp_fn(h)
                    local delta = torch.Tensor(data.beamDelta)
                    local energy = torch.Tensor(data.beamEnergyLimits)
                    local nel, np, temperature, HI, HII = data.nel, data.np, data.temperature, data.HI, data.HII
                    -- Loopwise position s is computed from xi and eta
                    -- local s = math.atan2(xi / d.aabb.xi.max, eta / d.aabb.eta.max) / (math.pi / 2)
                    -- local r = math.atan2(y / aabb.y.max, x / aabb.x.max) / (math.pi / 2)
                    local xx = math.floor(scale * x - scale * d.aabb.x.min + scale * d.trans.x)
                    local yy = math.floor(scale * y - scale * d.aabb.y.min + scale * d.trans.y)
                    local zz = math.floor(scale * z - scale * d.aabb.z.min + scale * d.trans.z)

                    local viewDir = math.pi/2

                    local b0x = cos(p.phi0) * p.b0Mag
                    local b0y = sin(p.phi0) * p.b0Mag

                    local xAxis = x - p.numVox / 2
                    local yAxis = y
                    local zAxis = z - p.numVox / 2

                    local b0 = maf.vec3(b0x, b0y, 0)
                    local rVec = maf.vec3(xAxis, yAxis, zAxis)
                    local mu = -b0 * p.height^3
                    local b = (rVec * (3 * mu:dot(rVec)) - mu * rVec:length()^2) / (rVec:length()^5)
                    local bMag = b:length()
                    local bDir = b:normalize():dot(d.rayDir)
                    local bAngle = math.acos(bDir)


                    if (x > d.aabb.x.min and x < d.aabb.x.max) and
                       (y > d.aabb.y.min and y < d.aabb.y.max) and
                       (z > d.aabb.y.min and z < d.aabb.z.max)
                    then
                        -- print(tostring(gyroIn))
                        -- print(("%x"):format(tonumber(ffi.cast('intptr_t', gyroIn))))
                        local idx = index(xx,yy,zz)
                        pool:addjob(
                            function() 
                                local gyroIn = C.GyroIn()
                                gyroIn.delta = delta:data()
                                gyroIn.deltaLen = delta:size(1)

                                gyroIn.energy = energy:data()
                                gyroIn.energyLen = energy:size(1)

                                gyroIn.frequency = frequency:data()
                                gyroIn.frequencyLen = frequency:size(1)

                                gyroIn.bMag = bMag
                                gyroIn.angle = bAngle * 180 / math.pi

                                if gyroIn.angle > 87 and gyroIn.angle < 90 then
                                    gyroIn.angle = 87
                                elseif gyroIn.angle >= 90 and gyroIn.angle < 93 then
                                    gyroIn.angle = 93
                                elseif gyroIn.angle > 177 and gyroIn.angle < 180 then
                                    gyroIn.angle = 177
                                elseif gyroIn.angle >= 180 and gyroIn.angle < 183 then
                                    gyroIn.angle = 183
                                elseif gyroIn.angle < 3 then
                                    gyroIn.angle = 3
                                end

                                gyroIn.nel = nel
                                gyroIn.np = np
                                local thermalData = C.ThermalRadiationData()
                                thermalData.temperature = temperature
                                thermalData.protonDensity = HII
                                thermalData.neutralHDensity = HI

                                local data = C.Gyro.GyroSimulateC(gyroIn, thermalData)
                                -- print(data.jo)
                                -- print(('%x'):format(tonumber(C.ffi.cast('intptr_t', data.jo))))
                                return torch.LongTensor({tonumber(C.ffi.cast('intptr_t', data.jo)),
                                            tonumber(C.ffi.cast('intptr_t', data.jx)),
                                            tonumber(C.ffi.cast('intptr_t', data.ko)),
                                            tonumber(C.ffi.cast('intptr_t', data.kx)),
                                            tonumber(C.ffi.cast('intptr_t', data.jtherm)),
                                            tonumber(C.ffi.cast('intptr_t', data.ktherm))})

                            end,
                            function(ret)
                                sampleCount[idx+1] = sampleCount[idx+1] + 1
                                local indices = {'jo', 'jx', 'ko', 'kx', 'jtherm', 'ktherm'}

                                if grid[idx].jo == nil then
                                    for i=1,#indices do 
                                        grid[idx][indices[i]] = ffi.cast('f64*', ret[i])
                                    end
                                else
                                    for i=1,#indices do 
                                        local coeff = ffi.cast('f64*', ret[i])
                                        for j = 0,frequency:size(1)-1 do
                                            grid[idx][indices[i]][j] = grid[idx][indices[i]][j] + coeff[j]
                                        end
                                    end
                                    -- ffi.C.free(ffi.cast('void*', ret[1]))
                                    freeList:append(ret[1])
                                end
                            end
                            )
                    end
                    arghDex = arghDex + 1
                elseif (x > d.aabb.x.min and x < d.aabb.x.max) and
                       (y > d.aabb.y.min and y < d.aabb.y.max) and
                       (z > d.aabb.y.min and z < d.aabb.z.max)
                then
                    local xx = math.floor(scale * x - scale * d.aabb.x.min + scale * d.trans.x)
                    local yy = math.floor(scale * y - scale * d.aabb.y.min + scale * d.trans.y)
                    local zz = math.floor(scale * z - scale * d.aabb.z.min + scale * d.trans.z)
                    local idx = index(xx,yy,zz)
                    emptySamples[idx+1] = emptySamples[idx+1] + 1
                end
            end
        end
        -- pool:synchronize()
    end
    pool:synchronize()

    for i = 1,#freeList do
        ffi.C.free(ffi.cast('void*', freeList[i]))
    end

    for x = 0, xRange-1 do
        for y = 0, yRange-1 do
            for z = 0, zRange-1 do
                local idx = index(x, y ,z)
                if grid[idx].jo ~= nil then
                    local indices = {'jo', 'jx', 'ko', 'kx', 'jtherm', 'ktherm'}
                    -- local sampleFactor = sampleCount:data()[idx] / (sampleCount:data()[idx] + emptySamples:data()[idx])
                    local sampleFactor = sampleCount[idx+1] + emptySamples[idx+1]
                    -- if emptySamples[idx] ~= 0 then sampleFactor = sampleFactor / emptySamples:data()[idx] end
                    for i = 1,#indices do
                        for f = 0,frequency:size(1)-1 do 
                            grid[idx][indices[i]][f] = grid[idx][indices[i]][f] / sampleFactor
                        end
                    end
                end
            end
        end
    end



    return {
            grid = grid,
            aabb = d.aabb,
            idx = index,
            scale = scale,
            params = p,
            rayDir = d.rayDir,
            trans = d.trans,
            rotMat = d.rotMat
        } 
end

--- Refine the Grid around the dipole's footpoints.
-- @tparam Grid d Grid returned from @{dipole}
-- @tparam GyroParamTable gyroParams the shared parameters for the Gyro calculation
-- @int scale the scale factor to apply
-- @number footpointFraction fraction of the y length of the AABB taken up by the footpoints
-- @tparam function(height)->PlasmaParams interp_fn function returning, for a given height, the values of N_el, N_p, temperature, N_HI, N_HII
-- @tparam threads.Threads Torch thread pool to use for computing the gyrosynchrotron coeffecients
-- @treturn List{Grid} a List containing the two refined Grids shrunk around the footpoints with their computed parameters inside
local function create_high_res_footpoints(d, gyroParams, scale, footpointFraction, interp_fn, pool)
    -- Create a first aabb around each footpoint, then refine
    local eps = 0.5
    local aabb = d.aabb
    local leftFootpointAabb = {
        x = {min = aabb.x.min + eps, max = 0.5 * (aabb.x.max + aabb.x.min) - eps},
        y = {min = aabb.y.min + eps, max = aabb.y.min + (aabb.y.max - aabb.y.min) * footpointFraction},
        z = {min = aabb.z.min + eps, max = aabb.z.max - eps}
    }

    local rightFootpointAabb = {
        x = {min = 0.5 * (aabb.x.max + aabb.x.min) + eps, max = aabb.x.max - eps},
        y = {min = aabb.y.min + eps, max = aabb.y.min + (aabb.y.max - aabb.y.min) * footpointFraction},
        z = {min = aabb.z.min + eps, max = aabb.z.max - eps}
    }

    local aabbs = {leftFootpointAabb, rightFootpointAabb}

    local floor, ceil = math.floor, math.ceil
    local min, max = math.min, math.max

    for i = 1,#aabbs do
        local a = aabbs[i]
        for x = a.x.min, a.x.max do
            for y = a.y.min, a.y.max do
                for z = a.z.min, a.z.max do
                    local xx = floor(x)
                    local yy = floor(y)
                    local zz = floor(z)
                    if d.grid[d.idx(xx, yy, zz)] ~= 0 then
                        a.x.min = min(a.x.min, xx)
                        a.x.max = max(a.x.max, xx)
                        a.y.min = min(a.y.min, yy)
                        a.y.max = max(a.y.max, yy)
                        a.z.min = min(a.z.min, zz)
                        a.z.max = max(a.z.max, zz)
                    end
                end
            end
        end
    end
    local a = leftFootpointAabb
    local eps2 = 1.5
    a.x.min = max(a.x.min - eps2, aabb.x.min + eps)
    a.x.max = min(a.x.max + eps2, 0.5 * (aabb.x.max + aabb.x.min) - eps) 
    a.y.min = max(a.y.min - eps2, aabb.y.min + eps)
    a.y.max = min(a.y.max + eps2, aabb.y.min + (aabb.y.max - aabb.y.min) * footpointFraction) 
    a.z.min = max(a.z.min - eps2, aabb.z.min + eps)
    a.z.max = min(a.z.max + eps2, aabb.z.max - eps) 
    a.xi = aabb.xi
    a.eta = aabb.eta

    a = rightFootpointAabb
    local eps2 = 1.5
    a.x.min = max(a.x.min - eps2, 0.5 * (aabb.x.max + aabb.x.min) + eps) 
    a.x.max = min(a.x.max + eps2, aabb.x.max - eps)
    a.y.min = max(a.y.min - eps2, aabb.y.min + eps)
    a.y.max = min(a.y.max + eps2, aabb.y.min + (aabb.y.max - aabb.y.min) * footpointFraction) 
    a.z.min = max(a.z.min - eps2, aabb.z.min + eps)
    a.z.max = min(a.z.max + eps2, aabb.z.max - eps) 
    a.xi = aabb.xi
    a.eta = aabb.eta
    -- aabbs should now be refined

    local leftParams = copy(d)
    leftParams.aabb = leftFootpointAabb
    -- local leftGrid = visualise_dipole_param(leftParams, scale, gyroParams, 'height', interp_fn)
    local leftGrid = dipole_in_aabb(leftParams, scale, gyroParams, interp_fn, pool)
    local rightParams = copy(d)
    rightParams.aabb = rightFootpointAabb
    -- local rightGrid = visualise_dipole_param(rightParams, scale, gyroParams, 'height', interp_fn)
    local rightGrid = dipole_in_aabb(rightParams, scale, gyroParams, interp_fn, pool)

    return List {leftGrid, rightGrid}
end

--- Returns the minimum of each element of two maf.vec3's.
local function min_per_elem(vec1, vec2)
    local min = math.min
    return maf.vec3(min(vec1.x, vec2.x),
                    min(vec1.y, vec2.y),
                    min(vec1.z, vec2.z)) 
end

--- Returns the maximum of each element of two maf.vec3's.
local function max_per_elem(vec1, vec2)
    local max = math.max
    return maf.vec3(max(vec1.x, vec2.x),
                    max(vec1.y, vec2.y),
                    max(vec1.z, vec2.z)) 
end

--- Returns the minimum element of a maf.vec3.
local function min_elem(vec)
    return math.min(vec.x, vec.y, vec.z)
end

--- Returns the maximum element of a maf.vec3.
local function max_elem(vec)
    return math.max(vec.x, vec.y, vec.z)
end

--- Compute the intersections of a ray and a box.
-- Based on Slab method presented in Siggraph 1998 by S. Owen. Added Tavian
-- Barnes' improvements using the IEE754 standard to handle edge cases that must
-- otherwise be handled explicitly. Minor improvements in the number of
-- comparisons with infinity.
-- @tparam AABB box the AABB enclosing the region of interest
-- @tparam maf.vec3 rayOrigin the origin of the ray
-- @tparam maf.vec3 rayDirection the normalised direction vector of the ray
-- @treturn number the first intersection t with the box such that rayOrigin + rayDirection * t = intersection
-- @treturn number the second intersection with the box
-- @treturn nil if the ray does not intersect with the box then nil is returned
local function box_ray_intersections(box, rayOrigin, rayDirection)
    local origin = rayOrigin 

    local minExtent = maf.vec3(box.x.min, box.y.min, box.z.min)
    local maxExtent = maf.vec3(box.x.max, box.y.max, box.z.max)

    local minOffset = minExtent - origin
    local maxOffset = maxExtent - origin

    local recipDir = maf.vec3(1.0 / rayDirection.x, 1.0 / rayDirection.y, 1.0 / rayDirection.z)

    local t1 = minOffset * recipDir
    local t2 = maxOffset * recipDir

    local tMin = min_per_elem(t1, t2) 
    local tMax = max_per_elem(t1, t2) 
    local inf = math.inf
    tMin = min_per_elem(tMin, maf.vec3(inf, inf, inf))
    tMax = max_per_elem(tMin, tMax)

    local min = max_elem(tMin)
    local max = min_elem(tMax)

    if max > min then
        return min, max
    end
end

--- Raymarch through the provided regions and return the average value for each ray.
-- @tparam Grid lowRes the basic unrefined grid valid everywhere in the dipole
-- @tparam List{Grid} highRes List of high resolution refined grids valid in certain locations (eg. footpoints). If none then provide an empty list.
-- @tparam int resolution the resolution of the image to trace
-- @tparam GyroParamTable gyroParams the shared parameters required for the Gyro calculation
-- @treturn DoubleCArray the array containing the output image
-- @treturn function(x,y)->index function used to index the array in x and y coordinates
local function average_path_trace(lowRes, highRes, resolution, gyroParams)
    local image = DoubleCArray((resolution+1)^2)

    local idx2 = uniform_2d_index_gen(resolution)
    local aux_sort = function (a, b) return a.min < b.min end

    local numVox = lowRes.params.numVox
    local PrimarySampleRate = 0.1
    local AuxSampleRate = 0.01
    local floor = math.floor
    local rayDir = lowRes.rayDir
    local rotMat= lowRes.rotMat
    for u = 1, resolution do
        for v = 1, resolution do
            local x = (u / resolution * numVox) - numVox / 2
            local y = (v / resolution * numVox) - numVox / 2
            local rayCount = 0

            local rayOrigin = rotMat * maf.vec3(x, y, 0.0) + maf.vec3(numVox / 2, numVox / 2, numVox / 2)
            local tMin, tMax = box_ray_intersections(lowRes.aabb, rayOrigin, rayDir)

            if tMin ~= nil and tMax ~= nil and tMin < tMax then
                local intersections = {primary = {min = tMin, max = tMax, region = lowRes}, auxiliary = List {}}
                for i = 1,#highRes do
                    local sMin, sMax = box_ray_intersections(highRes[i].aabb, rayOrigin, rayDir)
                    if sMin ~= nil and sMax ~= nil and sMin < sMax then
                        intersections.auxiliary:append({min = sMin, max = sMax, region = highRes[i]})
                    end
                end
                intersections.auxiliary:sort(aux_sort)

                local sortedIntList = List{}
                if #intersections.auxiliary > 0 then
                    local currentAux = 1
                    sortedIntList:append{min = tMin, max = intersections.auxiliary[currentAux].min, 
                                         region = intersections.primary.region, sampleRate = PrimarySampleRate}
                    while currentAux <= #intersections.auxiliary do
                        if currentAux ~= 1 then
                            sortedIntList:append{min = intersections.auxiliary[currentAux-1].max, 
                                                 max = intersections.auxiliary[currentAux].min, 
                                                 region = intersections.primary.region,
                                                 sampleRate = PrimarySampleRate}
                        end
                        sortedIntList:append(copy(intersections.auxiliary[currentAux]))
                        sortedIntList[#sortedIntList].sampleRate = AuxSampleRate
                        currentAux = currentAux + 1
                    end
                    sortedIntList:append({min = sortedIntList[#sortedIntList].max, max = tMax,
                                          region = intersections.primary.region, sampleRate = PrimarySampleRate})
                else
                    sortedIntList:append{min = tMin, max = tMax,
                                         region = intersections.primary.region, sampleRate = PrimarySampleRate}
                end

                local val = 0.0
                local lengthInPlasma = 0.0
                for j = #sortedIntList, 1, -1 do
                    local grid = sortedIntList[j].region

                    for t = sortedIntList[j].max - sortedIntList[j].sampleRate, 
                            sortedIntList[j].min, 
                            -sortedIntList[j].sampleRate 
                    do
                        local r = (rayOrigin + rayDir * t)
                        local length = sortedIntList[j].sampleRate

                        local xx = floor(grid.scale * (r.x - grid.aabb.x.min))
                        local yy = floor(grid.scale * (r.y - grid.aabb.y.min))
                        local zz = floor(grid.scale * (r.z - grid.aabb.z.min))

                        if (r.x > grid.aabb.x.min and r.x < grid.aabb.x.max) and
                           (r.y > grid.aabb.y.min and r.y < grid.aabb.y.max) and
                           (r.z > grid.aabb.y.min and r.z < grid.aabb.z.max)
                        then
                            if grid.grid[grid.idx(xx, yy, zz)] ~= 0.0 then
                                val = val + grid.grid[grid.idx(xx, yy, zz)] * length
                                lengthInPlasma = lengthInPlasma + length
                            end
                        end
                    end
                end
                local idx = idx2(u,v)
                if lengthInPlasma == 0.0 then
                    image[idx] = val
                else
                    image[idx] = val / lengthInPlasma
                end
            end
        end
    end
    return image, idx2
end

--- Raymarch through the provided regions and return the integrated brightness temperature values in O and X modes and thermal emission for each ray.
-- @tparam Grid lowRes the basic unrefined grid valid everywhere in the dipole
-- @tparam List{Grid} highRes List of high resolution refined grids valid in certain locations (eg. footpoints). If none then provide an empty list.
-- @tparam int resolution the resolution of the image to trace
-- @tparam GyroParamTable gyroParams the shared parameters required for the Gyro calculation
-- @treturn DoubleCArray the array containing the output image
-- @treturn function(x,y)->index function used to index the array in x and y coordinates
local function path_trace(lowRes, highRes, resolution, gyroParams)
    local imageList = List{}
    for i = 1, #gyroParams.frequency do
        local o = DoubleCArray((resolution+1)^2)
        local x = DoubleCArray((resolution+1)^2)
        local therm = DoubleCArray((resolution+1)^2)
        imageList:append({o = o, x = x, therm = therm})
    end

    local idx2 = uniform_2d_index_gen(resolution)
    local aux_sort = function (a, b) return a.min < b.min end

    local numVox = lowRes.params.numVox
    local PrimarySampleRate = 0.1
    local AuxSampleRate = 0.01
    local floor = math.floor
    local max = math.max
    local rayDir = lowRes.rayDir
    local rotMat = lowRes.rotMat
    for u = 1, resolution do
        for v = 1, resolution do
            local x = (u / resolution * numVox) - numVox / 2
            local y = (v / resolution * numVox) - numVox / 2
            local rayCount = 0

            local rayOrigin = rotMat * maf.vec3(x, y, 0.0) + maf.vec3(numVox / 2, numVox / 2, numVox / 2)
            local tMin, tMax = box_ray_intersections(lowRes.aabb, rayOrigin, rayDir)

            if tMin ~= nil and tMax ~= nil and tMin < tMax then
                local intersections = {primary = {min = tMin, max = tMax, region = lowRes}, auxiliary = List {}}
                for i = 1,#highRes do
                    local sMin, sMax = box_ray_intersections(highRes[i].aabb, rayOrigin, rayDir)
                    if sMin ~= nil and sMax ~= nil and sMin < sMax then
                        intersections.auxiliary:append({min = sMin, max = sMax, region = highRes[i]})
                    end
                end
                intersections.auxiliary:sort(aux_sort)

                local sortedIntList = List{}
                if #intersections.auxiliary > 0 then
                    local currentAux = 1
                    sortedIntList:append{min = tMin, max = intersections.auxiliary[currentAux].min, 
                                         region = intersections.primary.region, sampleRate = PrimarySampleRate}
                    while currentAux <= #intersections.auxiliary do
                        if currentAux ~= 1 then
                            sortedIntList:append{min = intersections.auxiliary[currentAux-1].max, 
                                                 max = intersections.auxiliary[currentAux].min, 
                                                 region = intersections.primary.region,
                                                 sampleRate = PrimarySampleRate}
                        end
                        sortedIntList:append(copy(intersections.auxiliary[currentAux]))
                        sortedIntList[#sortedIntList].sampleRate = AuxSampleRate
                        currentAux = currentAux + 1
                    end
                    sortedIntList:append({min = sortedIntList[#sortedIntList].max, max = tMax,
                                          region = intersections.primary.region, sampleRate = PrimarySampleRate})
                else
                    sortedIntList:append{min = tMin, max = tMax,
                                         region = intersections.primary.region, sampleRate = PrimarySampleRate}
                end

                -- local numSteps = 1

                local io = {} 
                for f = 1, #gyroParams.frequency do
                    io[#io + 1] = 0.0
                end
                local ix = copy(io)
                local itherm = copy(io)

                for j = #sortedIntList, 1, -1 do
                    local grid = sortedIntList[j].region

                    -- for t = sortedIntList[j].max - sortedIntList[j].sampleRate, 
                    for t = sortedIntList[j].max, 
                            sortedIntList[j].min, 
                            -sortedIntList[j].sampleRate 
                    do
                        local r = (rayOrigin + rayDir * t)
                        local length = sortedIntList[j].sampleRate * lowRes.params.VoxToCm

                        local xx = max(floor(grid.scale * (r.x - grid.aabb.x.min)), 0)
                        local yy = max(floor(grid.scale * (r.y - grid.aabb.y.min)), 0)
                        local zz = max(floor(grid.scale * (r.z - grid.aabb.z.min)), 0)

                        if (r.x > grid.aabb.x.min and r.x < grid.aabb.x.max) and
                           (r.y > grid.aabb.y.min and r.y < grid.aabb.y.max) and
                           (r.z > grid.aabb.y.min and r.z < grid.aabb.z.max)
                        then
                            -- print(xx, yy, zz)
                            if grid.grid[grid.idx(xx, yy, zz)].jo ~= nil then
                                local p = grid.grid[grid.idx(xx, yy, zz)]

                                for freq = 1, #gyroParams.frequency do
                                    local f = freq- 1;
                                    if p.ko[f] ~= 0.0 or p.ktherm[f] ~= 0.0 then
                                        io[freq] = io[freq] * math.exp(-(p.ko[f] + p.ktherm[f]) * length) 
                                                + p.jo[f] / (p.ko[f] + p.ktherm[f]) 
                                                  * (1 - math.exp(-(p.ko[f] + p.ktherm[f]) * length))
                                    else
                                        io[freq] = io[freq] + p.jo[f] * length
                                    end
                                    if p.kx[f] ~= 0.0 or p.ktherm[f] ~= 0.0 then
                                        ix[freq] = ix[freq] * math.exp(-(p.kx[f] + p.ktherm[f]) * length) 
                                                + p.jx[f] / (p.kx[f] + p.ktherm[f]) 
                                                  * (1 - math.exp(-(p.kx[f] + p.ktherm[f]) * length))
                                    else
                                        ix[freq] = ix[freq] + p.jx[f] * length
                                    end
                                    if p.ktherm[f] ~= 0.0 or p.ko[f] ~= 0.0 or p.kx[f] ~= 0.0 then
                                        itherm[freq] = itherm[freq] * math.exp(-(p.ko[f] + p.kx[f] + p.ktherm[f]) * length) 
                                                + p.jtherm[f] / (p.ko[f] + p.kx[f] + p.ktherm[f]) 
                                                  * (1 - math.exp(-(p.ko[f] + p.kx[f] + p.ktherm[f]) * length))
                                    else
                                        itherm[freq] = itherm[freq] + p.jtherm[f] * length
                                    end
                                end
                            end
                        end
                    end
                end
                -- print(u, v)
                -- grid2[idx2(u, v)] = (3e10)^2 * io / (2 * 1.38e-16 * gyroParams.frequency[frequencyIdx]^2)
                local idx = idx2(u,v)
                for f = 1, #gyroParams.frequency do
                    local ims = imageList[f]
                    local brightnessTempFactor = (3e10)^2 / (2 * 1.3806e-16 * gyroParams.frequency[f]^2)
                    ims.o[idx] = io[f] * brightnessTempFactor
                    ims.x[idx] = ix[f] * brightnessTempFactor
                    ims.therm[idx] = itherm[f] * brightnessTempFactor
                end
            end
        end
    end
    return imageList, idx2
end

--- Returns the rotation quaternion for a rotation about the z, then y, then x axes.
-- @number x the rotation to perform about the x axis in radians
-- @number y the rotation to perform about the y axis in radians
-- @number z the rotation to perform about the z axis in radians
-- @treturn maf.quat the quaternion defining this rotation
local function rotation_zyx(x, y, z)
    local rot = maf.rotation
    return rot():angleAxis(z, 0, 0, 1) * rot():angleAxis(y, 0, 1, 0) * rot():angleAxis(x, 1, 0, 0)
end

--- Returns the rotation quaternion for a given location on the Sun.
-- @number tilt the "lean" of the flare away from the normal to the surface in degrees
-- @number azimuth the rotation of the flare from the N-S axis in degrees
-- @number latitude the solar location in degrees
-- @number longitude the solar location in degrees
-- @treturn maf.quat the quaternion defining the rotation for this location
local function solar_location(tilt, azimuth, latitude, longitude)
    local til = tilt * math.pi / 180
    local az = azimuth * math.pi / 180
    local lat = latitude * math.pi / 180
    local lon = longitude * math.pi / 180

    return rotation_zyx(0, math.pi/2, 0) * rotation_zyx(math.pi/2, 0, 0) * rotation_zyx(0, -til, 0) * rotation_zyx(0, 0, az) * rotation_zyx(-lat, 0, 0) * rotation_zyx(0, -lon, 0) 
end

--- Returns the ray direction for a given rotation.
-- @tparam maf.quat rotMat the rotation to generate the direction vector for
-- @treturn maf.vec3 the normalised direction vector for this rotation
local function ray_direction(rotMat)
    return (rotMat * maf.vec3(0.0, 0.0, 1.0):normalize()):normalize()
end

--- Axis Aligned Bounding Box.
-- @table AABB
-- @tfield MinMax x
-- @tfield MinMax y
-- @tfield MinMax z

--- @table MinMax
-- @tfield number min
-- @tfield number max

--- @table exports 
-- @field ArcToCm conversion from arcseconds from Earth to cm on the surface of the Sun.
-- @export
return { solar_location = solar_location,
         rotation_zyx = rotation_zyx,
         ray_direction = ray_direction,
         path_trace = path_trace,
         average_path_trace = average_path_trace,
         visualise_dipole_param = visualise_dipole_param,
         create_high_res_footpoints = create_high_res_footpoints,
         dipole = dipole,
         create_gyro_thread_pool = create_gyro_thread_pool,
         dipole_in_aabb = dipole_in_aabb,
         backup_grid = backup_grid,
         restore_backup_grid = restore_backup_grid,
         DoubleCArray = DoubleCArray,
         ArcToCm = ArcToCm,
}