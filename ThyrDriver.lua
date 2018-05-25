--- Driver file for Thyr Gyrosynchrotron simulations
-- This file provides the specifics of the flare model we want to simulate, sets up and runs the simulation.
-- File to be read "C-style" i.e. with the main function at the bottom
assert(loadfile('Parse.lua'))()
local List = require('pl.List')
local Plyght = require('Plyght.Plyght')
local Thyr = require('Thyr')
local maf = require('maf')
local f1Data = parse_atmos_data('f1_adj.csv')
local c7Data = parse_atmos_data('c7_adj.csv')

local delta = {3}
local energyLimits = {10, 5000}

-- This is a simple method for finiding the top of the chromosphere as it just
-- looks at the provided atmosphere and takes 30 km below the maximum gradient
-- (in the transition region)
local function find_chromosphere_top(atmosData)
    local maxGrad = 0
    local maxIdx = 0
    for i = 2,#atmosData.height-1 do
        local grad = (atmosData.temperature[i-1] + atmosData.temperature[i+1]) / (atmosData.height[i+1] - atmosData.height[i-1])
        if grad > maxGrad then
            maxGrad = grad
            maxIdx = i
        end
    end
    return atmosData.height[maxIdx] - 3e6
end

-- Very primitive method for dropping off the non-thermal electron density as we
-- descend into the chromosphere. The non-thermal electron distribution (assumed
-- to be a power law) has a distribution of stopping depths. We scale the total
-- density as a function of height, rather than adjusting the lower bound of the
-- distribution as we descend into the atmosphere. Thus, there is room for
-- improvement in the way this is handled.
local function nel_stopping_array(atmosData)
    local maxChromoH = find_chromosphere_top(atmosData)
    local maxChromoNel = 1e8
    local coronaNel = 1e6
    local d = delta[1]
    local eMin = energyLimits[1]
    local eMax = energyLimits[2]

    atmosData.nel = {}
    for i = 1,#atmosData.height do
        if atmosData.height[i] > maxChromoH then
            atmosData.nel[i] = coronaNel
        else
            atmosData.nel[i] = maxChromoNel
        end
    end

    local function create_log_space(eMin, eMax, numPts)
        local lMin = math.log10(eMin)
        local lMax = math.log10(eMax)
        local lStep = (lMax - lMin) / (numPts - 1)
        local space = {}
        for j = 1,numPts do
            space[j] = 10^(lMin + lStep * (j-1))
        end
        return space
    end

    -- Create the power law probabilities
    local function ElectronTable() return {energy = 0, prob = 0, columnDensity = 0, stopHeight = 0} end
    local lSpace = create_log_space(eMin, eMax, 100)
    local electrons = {}
    for j = 1,#lSpace do
        electrons[j] = ElectronTable()
        electrons[j].energy = lSpace[j]
        electrons[j].prob = lSpace[j]^(-d)
    end
    -- Normalise
    local sum = 0.0
    for _,v in ipairs(electrons) do
        sum = sum + v.prob
    end
    for _,v in ipairs(electrons) do
        v.prob = v.prob / sum
    end
    -- Find the top of the chromosphere where our max density is
    local topChromoI = 1
    while topChromoI < #atmosData.height do
        if atmosData.height[topChromoI] > maxChromoH then
            break
        end
        topChromoI = topChromoI + 1
    end

    -- Get the average density for step and compute the total column density for the electron
    for _,pt in ipairs(electrons) do
        for j = topChromoI, 2, -1 do
            -- local density = atmosData.np[j]
            -- local densityP1 = atmosData.np[j-1]

            local density = atmosData.HI[j] + atmosData.HII[j]
            local densityP1 = atmosData.HI[j-1] + atmosData.HII[j-1]
            local aveDensity = 0.5 * (density + densityP1)
            local deltaH = math.abs(atmosData.height[j-1] - atmosData.height[j])
            pt.columnDensity = pt.columnDensity + aveDensity * deltaH
            -- If the column density is above the stop density then input the stop height and move onto the next step
            if pt.columnDensity >= 1e17 * pt.energy * pt.energy then
                pt.stopHeight = atmosData.height[j-1]
                break
            end
        end
    end

    -- Reduce the electron density
    for _,pt in ipairs(electrons) do
        for hIndex,h in ipairs(atmosData.height) do
            if h < pt.stopHeight then
                atmosData.nel[hIndex] = atmosData.nel[hIndex] - pt.prob * maxChromoNel
            end
        end
    end
end

local function plot_atmos_data(atmosData, prefix)
    prefix = prefix and prefix ..'.png' or 'atmosData.png'
    Plyght:start_frame()
          :plot()
          :x_range(0, 5e8)
          :plot_type('semilogy')
          :line_label('Temperature')
          :x_label('Height [cm]')
          :y_label('Temperature [K]')
          :line(atmosData.height, atmosData.temperature)
          :legend()
          :plot()
          :x_range(0, 5e8)
          :plot_type('semilogy')
          :x_label('Height [cm]')
          :y_label('Number density [cm$^-3$]')
          :line_label('Non-thermal electron density')
          :line(atmosData.height, atmosData.nel)
          :line_label('Neutral hydrogen density')
          :line(atmosData.height, atmosData.HI)
          :line_label('Proton density')
          :line(atmosData.height, atmosData.HII)
          :legend()
          :print(prefix)
          :end_frame()
end

function main()
    -- Initialise
    local startTime = os.time()
    Plyght:init()
    -- Load in atmospheric model
    local atmosData = c7Data
    local prefix = 'c7Data'
    -- Compute non-thermal electron drop off with increasing density
    nel_stopping_array(atmosData)
    -- Produce a plot of the atmosphere if Plyght is running
    plot_atmos_data(atmosData, prefix)

    -- Number of low-resolution voxels to be used per side of cube.
    local numVox = 32
    -- Side length of voxel grid, i.e. this needs to be able to contain the flare.
    local gridSide = 80 * Thyr.ArcToCm
    local VoxToCm = gridSide / numVox
    local height = numVox / 1.5
    local dipoleParams = 
    {
        numVox = numVox,
        height = height,
        depth = numVox / 5,
        phi0 = 0,
        rho0 = height / 5,
        gridSide = gridSide,
        b0Mag = 100,
        VoxToCm = VoxToCm,
    }
    print(("Voxel in km: %f, \": %f"):format(VoxToCm / 1e5, VoxToCm / Thyr.ArcToCm ))
    print(("Height in km: %f, \": %f"):format(height * VoxToCm / 1e5, height * VoxToCm / Thyr.ArcToCm ))

    -- Make basic dipole shape
    local resolution = 512 
    local trans = maf.vec3(0, 0, 0)
    print('Making Basic Dipole: ', os.time() - startTime)
    local grid3, idx, aabb, count = Thyr.dipole(dipoleParams, trans)
    dipoleParams.count = count

    -- Scale factor for high res footpoints (scales voxel side length)
    local scale  = 1
    -- Create interpolation function for the atmospheric data
    local interp_fn = function(height) return interpolate_data(height, atmosData) end
    print('Making Accurate Dipole: ', os.time() - startTime)
    -- Create the rotation quaternion for the object (effectively just a fancy form of a rotation matrix)
    local rotMat = Thyr.solar_location(0, -20, 30, 70)
    -- Deduce the unit direction vector of our parallel rays from the quaternion
    local rayDir = Thyr.ray_direction(rotMat)

    -- Parameters describing the dipole
    local dipoleData = {grid = grid3, 
                        aabb = aabb, 
                        idx = idx, 
                        scale = 1, 
                        params = dipoleParams, 
                        rayDir = rayDir, 
                        trans = trans,
                        rotMat = rotMat,
                    }
    -- Parameters needed for the simulation of GS emission
    local gyroParams = { frequency = {34e9}, energy = energyLimits, delta = delta }
    local plot

    if true then
        -- Perform gyro simulation
        local grid3HD = Thyr.dipole_in_aabb(dipoleData, scale, gyroParams, interp_fn)
        print('Making High-Res Footpoints: ', os.time() - startTime)
        local highResRegions = Thyr.create_high_res_footpoints(dipoleData, gyroParams, 1, 0.25, interp_fn)

        if false then
            -- This code shows how to serialise the simulation data to disk and load it again
            Thyr.backup_grid(grid3HD, gyroParams, prefix..'/base')
            Thyr.backup_grid(highResRegions[1], gyroParams, prefix..'/hr1')
            Thyr.backup_grid(highResRegions[2], gyroParams, prefix..'/hr2')

            local grid3HD2, gyroParams2 = Thyr.restore_backup_grid(prefix..'/base')
            local hr1, _ = Thyr.restore_backup_grid(prefix..'/hr1')
            local hr2, _ = Thyr.restore_backup_grid(prefix..'/hr2')
            local highRes = List {hr1, hr2}
        end

        print('Path Tracing: ', os.time() - startTime)
        local imageList, idx2 = Thyr.path_trace(grid3HD, highResRegions, resolution, gyroParams)

        plot = function(mode, freqIdx)
            Plyght:start_frame()
                :plot()
                :colorbar()
                :imshow(imageList[freqIdx][mode], resolution, resolution)
                :end_frame()
        end
        plot('o', 1)
    else
        -- Display one of the parameters of the simulations parameters averaged along each line of sight
        local grid3HD = Thyr.visualise_dipole_param(dipoleData, scale, gyroParams, 'height', interp_fn)
        local highRes = List{}
        local image, idx2 = Thyr.average_path_trace(grid3HD, highRes, resolution, gyroParams)

        print('Plotting: ', os.time() - startTime)

        plot = function(mode, freqIdx)
            Plyght:start_frame()
                :plot()
                :colorbar()
                :imshow(image, resolution, resolution)
                :end_frame()
        end
        plot()
    end

    -- return the plotting function so it can be used interactively, in the case
    -- where it takes different modes/frequencies to visualise the output
    return plot
    -- Don't close Plyght if you want to use this interactively
    -- Plyght:close()
end

plot = main()