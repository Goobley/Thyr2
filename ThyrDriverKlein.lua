assert(loadfile('Parse.lua'))()
local List = require('pl.List')
local stringx = require('pl.stringx')
stringx.import()
local Thyr = require('Thyr')
local maf = require('maf')
local Plyght = require('Plyght.Plyght')
local f1Data = parse_atmos_data('f1_adj.csv')
local c7Data = parse_atmos_data('c7_adj.csv')

local delta = {3}
local energyLimits = {30, 1000}

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
    return atmosData.height[maxIdx] - 5e6
end

local function nel_stopping_array(atmosData)
    -- Drop off non-thermal electrons with stopping column density as we descend
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

    -- Plyght:start_frame()
    --       :plot()
    --       :plot_type('semilogy')
    --       :x_range(0, 1.5 * maxChromoH)
    --       :line(atmosData.height, atmosData.nel)
    --       :end_frame()

end

local function plot_atmos_data(atmosData, prefix)
    prefix = prefix and prefix ..'.png' or 'atmosData.png'
    local h = find_chromosphere_top(atmosData)
    Plyght:start_frame()
          :plot()
          :x_range(0, 10 * h)
          :plot_type('semilogy')
          :line_label('Temperature')
          :x_label('Height [cm]')
          :y_label('Temperature [K]')
          :line(atmosData.height, atmosData.temperature)
          :legend()
          :plot()
          :x_range(0, 10 * h)
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
    local startTime = os.time()
    Plyght:init()
    local atmosData = f1Data
    -- local prefix = 'KleinM5pc'
    -- local prefix = 'Klein30'
    local prefix = 'Klein50P5pc'
    nel_stopping_array(atmosData)
    plot_atmos_data(atmosData, prefix)

    local numVox = 256
    -- local gridSide = 80 * Thyr.ArcToCm
    -- local height = numVox / 1.5
    local gridSide = (2e10 * 1.5)  * 2
    local VoxToCm = gridSide / numVox
    local height = 1.9e10 / VoxToCm + (numVox / 3)
    local dipoleParams = 
    {
        numVox = numVox,
        height = height,
        depth = numVox / 2.86,
        -- phi0 = math.pi / 8,
        phi0 = 0,
        -- rho0 = height / 5,
        rho0 = 0.1e10 / VoxToCm,
        gridSide = gridSide,
        b0Mag = 3.21 * 1.05,
        -- b0Mag = 3.3,
        VoxToCm = VoxToCm,
    }
    print(("Voxel in km: %f, \": %f"):format(VoxToCm / 1e5, VoxToCm / Thyr.ArcToCm ))
    print(("Height in km: %f, \": %f"):format(height * VoxToCm / 1e5, height * VoxToCm / Thyr.ArcToCm ))

    local resolution = 512
    local trans = maf.vec3(0, 0, 0)
    print('Making Basic Dipole: ', os.time() - startTime)
    -- if true then
        local grid3, idx, aabb, count = Thyr.dipole(dipoleParams, trans)
        dipoleParams.count = count

        local scale  =  1
        -- local interp_fn = function(height) return interpolate_data(height, atmosData) end
        -- The Klein loop is large compared to the curvature of the Sun, so the
        -- footpoints end up at a noticeably negative altitude compared to the
        -- altitude of the surface. For sensible loops like we have the ability
        -- to observe in Solar Physics these days, this locally planar
        -- assumption is fine. Here however, it isn't. So let's apply a fudge
        -- factor to compensate.
        local interp_fn = function(height) 
            height = height + 1.52e9
            local n0 = 1e7
            -- local n = n0 * math.exp((1.95e10 - height) / (((height / Thyr.ArcToCm + 960) / 960)^2 * 1.19e10))
            local n = n0 * math.exp((2e10 - height) / (((height / Thyr.ArcToCm + 960) / 960)^2 * 1.02e10))
            return {height = height,
                    temperature = 1,
                    np = n,
                    HI = 1,
                    HII = 1,
                    HMinus = 1,
                    nel = 1e4,
                    beamDelta = delta,
                    beamEnergyLimits = energyLimits }
        end
        print('Making Accurate Dipole: ', os.time() - startTime)
        -- local rotMat = rotation_zyx(math.pi/4, 0, 0):mul(rotation_zyx(math.pi/2, 0, 0):mul(rotation_zyx(0, 0, math.pi/2)))
        -- local rotMat = Thyr.solar_location(0, -20, 30, 70)
        local rotMat
        if prefix:startswith('Klein50') then
            rotMat = Thyr.solar_location(0, 90, 0, 50)
        elseif prefix:startswith('Klein30') then
            rotMat = Thyr.solar_location(0, 90, 0, 30)
        else
            rotMat = Thyr.solar_location(0, 90, 0, 0)
        end
        -- rotMat = Thyr.solar_location(90, 0, 0, 0)
        local rayDir = Thyr.ray_direction(rotMat)

        local dipoleData = {grid = grid3, 
                            aabb = aabb, 
                            idx = idx, 
                            scale = 1, 
                            params = dipoleParams, 
                            rayDir = rayDir, 
                            trans = trans,
                            rotMat = rotMat,
                        }
        local gyroParams = { 
            frequency = {80e6, 160e6, 320e6},
                                          energy = energyLimits, delta = delta }
        -- local grid3HD = Thyr.visualise_dipole_param(dipoleData, scale, gyroParams, 'N_p', interp_fn)

        local gyroThreads = 4
        local gyroFile = io.open('CoreCount', 'r')
        if gyroFile ~= nil then gyroThreads = gyroFile:read('*n') gyroFile:close() end


        local pool = Thyr.create_gyro_thread_pool(gyroThreads)
        local grid3HD = Thyr.dipole_in_aabb(dipoleData, scale, gyroParams, interp_fn, pool)
        -- print('Making High-Res Footpoints: ', os.time() - startTime)
        -- local highResRegions = Thyr.create_high_res_footpoints(dipoleData, gyroParams, 8, 0.25, interp_fn)
    -- end


    Thyr.backup_grid(grid3HD, gyroParams, prefix..'/base')
    -- -- Thyr.backup_grid(highResRegions[1], gyroParams, prefix..'/hr1')
    -- -- Thyr.backup_grid(highResRegions[2], gyroParams, prefix..'/hr2')

    local grid3HD2, gyroParams2 = Thyr.restore_backup_grid(prefix..'/base')
    -- -- local hr1, _ = Thyr.restore_backup_grid(prefix..'/hr1')
    -- -- local hr2, _ = Thyr.restore_backup_grid(prefix..'/hr2')
    -- -- local highRes = List {hr1, hr2}
    local highRes = List{}
    -- -- local highRes = highResRegions

    print('Path Tracing: ', os.time() - startTime)

    local imageList, idx2 = Thyr.path_trace(grid3HD2, highRes, resolution, gyroParams2)

    print('Plotting: ', os.time() - startTime)

    local plot = function(mode, freqIdx)
        Plyght:start_frame()
            :plot()
            :colorbar()
            :imshow(imageList[freqIdx][mode], resolution, resolution)
            :end_frame()
    end

    -- local image, idx2 = Thyr.average_path_trace(grid3HD, highRes, resolution, gyroParams2)

    -- print('Plotting: ', os.time() - startTime)

    -- local plot = function(mode, freqIdx)
    --     Plyght:start_frame()
    --         :plot()
    --         :colorbar()
    --         :imshow(image, resolution, resolution)
    --         :end_frame()
    -- end
    -- plot()

    -- return plot
    -- Plyght:close()

end

plot = main()