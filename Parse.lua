function LerpHeight(height, heightRange)
   local xs = heightRange
   return function(ys)
      return ys[1] + (height - xs[1]) * (ys[2] - ys[1]) / (xs[2] - xs[1])
   end
end

function ParseCSVLine(line)
   -- Handles numbers only!!
   local res = {}
   local pos = 1
   sep = ','
   while true do
      local char = string.sub(line, pos, pos)
      local startPt, endPt = string.find(line, sep, pos)
      if (startPt) then
         table.insert(res, tonumber(string.sub(line, pos, startPt-1)))
         pos = endPt + 1
      else
         table.insert(res, tonumber(string.sub(line, pos)))
         break
      end
   end
   return res
end

function parse_atmos_data(filename)
   -- init data
   local data = {}
   data.height = {}
   data.temperature = {}
   data.np = {}
   data.HI = {}
   data.HII = {}
   data.HMinus = {}

   -- read file and file table
   local f = assert(io.open(filename, 'r'))
   -- ignore first line with headers
   f:read('*line')
   while true do
      local line = f:read('*line')
      if line == nil then break end
      local csv = ParseCSVLine(line)
      table.insert(data.height, csv[1])
      table.insert(data.temperature, csv[2])
      table.insert(data.np, csv[3])
      table.insert(data.HI, csv[4])
      table.insert(data.HII, csv[5])
      table.insert(data.HMinus, csv[6])
   end
   f:close()

   -- TODO(Chris): Maybe tidy up later by adding a method that returns the table on its own, for now we just copy and paste
   -- data.LerpExtract = function()

   return data
end

function parse_two_col_csv(filename)
    local xs, ys = {}, {}

    local f = assert(io.open(filename, 'r'))
    while true do
        local line = f:read('*line')
        if line == nil then break end
        local csv = ParseCSVLine(line)
        xs[#xs+1] = csv[1]
        ys[#ys+1] = csv[2]
    end
    f:close()
    return xs, ys
end

function FindClosestData(height, data)
   local h = data.height;
   for i = 2, #data.height do
      if height < h[i] then
         if math.abs(height - h[i]) < math.abs(height - h[i-1]) then
            return {height = data.height[i],
                    temperature = data.temperature[i],
                    np = data.np[i],
                    HI = data.HI[i],
                    HII = data.HII[i],
                    HMinus = data.HMinus[i]}
         else
            return {height = data.height[i-1],
                    temperature = data.temperature[i-1],
                    np = data.np[i-1],
                    HI = data.HI[i-1],
                    HII = data.HII[i-1],
                    HMinus = data.HMinus[i-1]}
         end
      end
   end
   return nil;
end

function interpolate_data(height, data)
   local h = data.height
   -- Check if factor is in array
   for i = 2, #data.height do
      if height > h[i-1] and height <= h[i] then
         local Lerp = LerpHeight(height, {data.height[i-1], data.height[i]})

         return {height = height,
                 temperature = Lerp({data.temperature[i-1], data.temperature[i]}),
                 np = Lerp({data.np[i-1], data.np[i]}),
                 HI = Lerp({data.HI[i-1], data.HI[i]}),
                 HII = Lerp({data.HII[i-1], data.HII[i]}),
                 HMinus = Lerp({data.HMinus[i-1], data.HMinus[i]})}
      end
   end

   -- If it isn't then use an end value
   if height > h[#h] then
      return {height = height,
              temperature = data.temperature[#h],
              np = data.np[#h],
              HI = data.HI[#h],
              HII = data.HII[#h],
              HMinus = data.HMinus[#h]}
   else
      return {height = height,
              temperature = data.temperature[1],
              np = data.np[1],
              HI = data.HI[1],
              HII = data.HII[1],
              HMinus = data.HMinus[1]}
   end
   -- in case something goes wrong
   return nil
end
