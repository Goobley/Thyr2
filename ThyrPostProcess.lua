assert(loadfile('Parse.lua'))()
local Thyr = require('Thyr')
local Plyght = require('Plyght.Plyght')
local List = require('pl.List')


local CgsToSfu = 1.0e19
local AstroUnit = 1.49597870e13
local ArcToCMSun = math.pi / 180.0 / 3600.0 * AstroUnit
local SpeedLight = 3e10
local Boltzmann = 1.3806e-16

local Dpi = 300

local function compute_total_emission_maps(imList, idx, frequencies, resolution)
    local totalMaps = List()
    for f = 1, #frequencies do
        totalMaps:append(Thyr.DoubleCArray((resolution+1)^2))
        for u = 1, resolution do
            for v = 1, resolution do
                local i = idx(u, v)
                totalMaps[f][i] = imList[f]['o'][i] + imList[f]['x'][i] + imList[f]['therm'][i]
            end
        end
    end
    return totalMaps
end

local function integrate_flux(imList, idx, resolution, viewSize, frequencies)
    local pxArea = (viewSize / resolution)^2
    local pxSr = pxArea / AstroUnit^2
    local flux = List()
    for f = 1, #frequencies do
        local intens = 0
        for u = 1, resolution do
            for v = 1, resolution do
                local i = idx(u, v)
                intens = intens + imList[f]['o'][i] + imList[f]['x'][i] + imList[f]['therm'][i]
            end
        end

        local tbToI = 2 * Boltzmann * frequencies[f]^2 / SpeedLight^2
        local fluxFactor = tbToI * pxSr * CgsToSfu
        intens = intens * fluxFactor
        flux:append(intens)
    end
    return flux
end

local function integrate_flux_mode(mode, imList, idx, resolution, viewSize, frequencies)
    local pxArea = (viewSize / resolution)^2
    local pxSr = pxArea / AstroUnit^2
    local flux = List()
    for f = 1, #frequencies do
        local intens = 0
        for u = 1, resolution do
            for v = 1, resolution do
                local i = idx(u, v)
                intens = intens + imList[f][mode][i]
            end
        end

        local tbToI = 2 * Boltzmann * frequencies[f]^2 / SpeedLight^2
        local fluxFactor = tbToI * pxSr * CgsToSfu
        intens = intens * fluxFactor
        flux:append(intens)
    end
    return flux
end

local function integrate_bg_flux(imList, idx, resolution, viewSize, frequencies)
    local pxArea = (viewSize / resolution)^2
    local pxSr = pxArea / AstroUnit^2
    local flux = List()
    for f = 1, #frequencies do
        local intens = 0
        for u = 1, resolution do
            for v = 1, resolution do
                local i = idx(u, v)
                intens = imList[f]['therm'][i]
            end
        end

        local tbToI = 2 * Boltzmann * frequencies[f]^2 / SpeedLight^2
        local fluxFactor = tbToI * pxSr * CgsToSfu
        intens = intens * fluxFactor
        flux:append(intens)
    end
    return flux
end

local function b_field_angle(grid, gyroParams, prefix)
    local atmosData
    if prefix == 'f1Data' then
        atmosData = parse_atmos_data('f1_adj.csv')
    else
        atmosData = parse_atmos_data('c7_adj.csv')
    end
    atmosData.nel = {}
    for i = 1, #atmosData.height do
        atmosData.nel[i] = 0.0
    end
    local interp_fn = function(height) return interpolate_data(height, atmosData) end
    local bAngle = Thyr.visualise_dipole_param(grid, 1, gyroParams, 'angle', interp_fn)
    return bAngle
end

local function load_grids_and_trace(prefix, resolution)
    local grid3HD2, gyroParams2 = Thyr.restore_backup_grid(prefix..'/base')
    local hr1, _ = Thyr.restore_backup_grid(prefix..'/hr1')
    local hr2, _ = Thyr.restore_backup_grid(prefix..'/hr2')
    local highRes = List {hr1, hr2}

    local imageList, idx2 = Thyr.path_trace(grid3HD2, highRes, resolution, gyroParams2)
    return imageList, idx2, grid3HD2, gyroParams2, highRes
end

local function signum(x)
    if x >= 0 then return 1
    else return -1 end
end

local function compute_polarisation_fraction_maps(imList, idx, bAngle, frequencies, resolution)
    local polMaps = List()
    local cos = math.cos
    local DToR = math.pi / 180
    for f = 1, #frequencies do
        polMaps:append(Thyr.DoubleCArray((resolution+1)^2))
        for u = 1, resolution do
            for v = 1, resolution do
                local i = idx(u, v)
                local angleFact = signum(cos(bAngle[i] * DToR))
                local totEmission = imList[f]['o'][i] + imList[f]['x'][i] + imList[f]['therm'][i]
                if totEmission == 0 then
                    polMaps[f][i] = 0
                else
                    polMaps[f][i] = angleFact * (imList[f]['o'][i] - imList[f]['x'][i]) / totEmission
                end
            end
        end
    end
    return polMaps
end

local function map_to_csv(prefix, imList, idx, frequencies, resolution)
    for f = 1, #frequencies do
        local csv = io.open(prefix..'_'..f..'.csv', 'w') 
        for u = 1, resolution do
            for v = 1, resolution do
                local i = idx(u, v)
                csv:write(imList[f][i])
                if v ~= resolution then
                    csv:write(',')
                end
            end
            if u ~= resolution then
                csv:write('\n')
            end
        end
        csv:flush()
        csv:close()
    end
end

local function totemission_to_csv(prefix, imList, idx, frequencies, resolution)
    for f = 1, #frequencies do
        local csv = io.open(prefix..'_tot_'..f..'.csv', 'w') 
        for u = 1, resolution do
            for v = 1, resolution do
                local i = idx(u, v)
                csv:write(imList[f]['o'][i] + imList[f]['x'][i] + imList[f]['therm'][i])
                if v ~= resolution then
                    csv:write(',')
                end
            end
            if u ~= resolution then
                csv:write('\n')
            end
        end
        csv:flush()
        csv:close()
    end
end

function main()
    Plyght:init()
    Plyght:start_frame():plot():fig_size(6, 6):end_frame()
    local prefix = 'c7DataHR'
    local title = 'High-Res C7 Atmosphere'
    local bgPrefix = 'c7Data'

    local resolution = 512

    local startTime = os.time()
    print('Restoring grids and path tracing: ', os.time() - startTime)
    local imageList, idx2, grid3HD2, gyroParams2, highRes = load_grids_and_trace(prefix, resolution)

    print('Summing Maps: ', os.time() - startTime)

    local totalMaps = compute_total_emission_maps(imageList, idx2, gyroParams2.frequency, resolution)

    print('Plotting Brightness Temperature: ', os.time() - startTime)

    local plot = function(freqIdx)
        Plyght:start_frame()
            :plot()
            :colorbar()
            :fig_size(6,6)
            :title(('Total Brightness Temperature at %.2f GHz!!nfor '):format(gyroParams2.frequency[freqIdx] / 1e9) .. title)
            :imshow(totalMaps[freqIdx], resolution, resolution)
            :print('TotTb_'..prefix..'_'..freqIdx..'.png', Dpi)
            :end_frame()
    end
    local plot2 = function()
        for i = 1, #totalMaps do
            plot(i)
        end
    end
    plot2()

    print('Plotting Integrated Flux: ', os.time() - startTime)
    local integratedFlux = integrate_flux(imageList, idx2, resolution, grid3HD2.params.VoxToCm * grid3HD2.params.numVox, gyroParams2.frequency)
    local oFlux = integrate_flux_mode('o', imageList, idx2, resolution, grid3HD2.params.VoxToCm * grid3HD2.params.numVox, gyroParams2.frequency)
    local xFlux = integrate_flux_mode('x', imageList, idx2, resolution, grid3HD2.params.VoxToCm * grid3HD2.params.numVox, gyroParams2.frequency)
    local thermFlux = integrate_flux_mode('therm', imageList, idx2, resolution, grid3HD2.params.VoxToCm * grid3HD2.params.numVox, gyroParams2.frequency)
    local maxIntFlux = 0
    for k = 1,#integratedFlux do
        maxIntFlux = math.max(maxIntFlux, integratedFlux[k]);
    end

    local plot3 = function()
        Plyght:start_frame()
            :plot()
            :fig_size(6,4)
            :x_label('Frequency [Hz]')
            :y_label('Integrated Flux [sfu]')
            :title('Integrated Flux vs. Frequency for '..title)
            :plot_type('loglog')
            :line_label('Total Flux')
            :line(gyroParams2.frequency, integratedFlux)
            :line_label('X-mode Flux')
            :line(gyroParams2.frequency, xFlux)
            :line_label('O-mode Flux')
            :line(gyroParams2.frequency, oFlux)
            :line_label('Thermal Flux')
            :line(gyroParams2.frequency, thermFlux)
            :y_range(1, 1.1*maxIntFlux)
            :legend()
            :print('IntFlux_'..prefix..'.png', Dpi)
            :end_frame()
    end
    plot3()

    if prefix ~= bgPrefix then
        print('Plotting Integrated Excess: ', os.time() - startTime)
        local bgIm, _, bgGrid, bgParams, _ = load_grids_and_trace(bgPrefix, resolution)
        local bgFlux = integrate_bg_flux(bgIm, idx2, resolution, bgGrid.params.VoxToCm * bgGrid.params.numVox, bgParams.frequency)
        local excessFlux = List()
        for i = 1, #bgParams.frequency do
            excessFlux:append(integratedFlux[i] - bgFlux[i])
        end
        Plyght:start_frame()
            :plot()
            :fig_size(6,4)
            :x_label('Frequency [Hz]')
            :y_label('Integrated Flare Excess [sfu]')
            :title('Integrated Excess vs. Frequency for '..title)
            :plot_type('loglog')
            :line(gyroParams2.frequency, excessFlux)
            :print('IntExcess_'..prefix..'.png', Dpi)
            :end_frame()
    end


    print('Plotting Polarisation Fraction: ', os.time() - startTime)
    local bAngle = b_field_angle(grid3HD2, gyroParams2, prefix)
    local bAngleImage = Thyr.average_path_trace(bAngle, List(), resolution, gyroParams2)

    local polMaps = compute_polarisation_fraction_maps(imageList, idx2, bAngleImage, gyroParams2.frequency, resolution)
    local plot4 = function(freqIdx)
        local range = 0
        local abs, max = math.abs, math.max
        for u = 1, resolution do
            for v = 1, resolution do
                local i = idx2(u,v)
                range = max(abs(polMaps[freqIdx][i]), range)
            end
        end
        Plyght:start_frame()
            :plot()
            :fig_size(6,6)
            :colorbar()
            :colormap('seismic')
            :title(('Polarisation Fraction at %.2f GHz!!nfor '):format(gyroParams2.frequency[freqIdx] / 1e9) .. title)
            :imshow(polMaps[freqIdx], resolution, resolution, -range, range)
            :print('PolFrac_'..prefix..'_'..freqIdx..'.png', Dpi)
            :end_frame()
    end
    local plotPol = function()
        for i = 1, #polMaps do
            plot4(i)
            -- print(gyroParams2.frequency[i])
            -- io.read()
        end
    end
    plotPol()

    totemission_to_csv(prefix, imageList, idx2, gyroParams2.frequency, resolution)
    map_to_csv(prefix..'_pol', polMaps, idx2, gyroParams2.frequency, resolution)

    return plotPol


end

plot = main()