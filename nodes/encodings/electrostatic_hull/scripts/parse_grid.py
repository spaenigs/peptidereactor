#!/usr/bin/env python

__author__ = "Cedric Staniewski"

#
# File Format
#
# OpenDX scalar data
# http://www.poissonboltzmann.org/docs/file-format-info/#opendx-scalar-data
#
#   object 1 class gridpositions counts nx ny nz
#   origin xmin ymin zmin
#   delta hx 0.0 0.0
#   delta 0.0 hy 0.0
#   delta 0.0 0.0 hz
#   object 2 class gridconnections counts nx ny nz
#   object 3 class array type double rank 0 items n data follows
#   u(0,0,0) u(0,0,1) u(0,0,2)
#   ...
#   u(0,0,nz-3) u(0,0,nz-2) u(0,0,nz-1)
#   u(0,1,0) u(0,1,1) u(0,1,2)
#   ...
#   u(0,1,nz-3) u(0,1,nz-2) u(0,1,nz-1)
#   ...
#   u(0,ny-1,nz-3) u(0,ny-1,nz-2) u(0,ny-1,nz-1)
#   u(1,0,0) u(1,0,1) u(1,0,2)
#   ...
#   attribute "dep" string "positions"
#   object "regular positions regular connections" class field
#   component "positions" value 1
#   component "connections" value 2
#   component "data" value 3

# # for python2 compatibility
# from __future__ import print_function

from collections import namedtuple
import copy, re, sys

file_header = []
file_footer = []

# named tuples
OpenDX = namedtuple('OpenDX', ['data', 'origin', 'delta', 'file_header', 'file_footer'])
class Point(namedtuple('Point', ['x', 'y', 'z'])):
  __slots__ = ()

  @property
  def label(self):
    return ",".join([str(x) for x in [self.x, self.y, self.z]])

  def add(self, p):
    return Point(self.x + p.x, self.y + p.y, self.z + p.z)

  def multiply(self, p):
    return Point(self.x * p.x, self.y * p.y, self.z * p.z)


def writeDX(dx, filename=''):
  f = sys.stdout
  if filename:
    f = open(filename, 'w')

  try:
    print(*dx.file_header, sep='', end='', file=f)
    for x in range(len(dx.data)):
      for y in range(len(dx.data[x])):
        for z in range(0, len(dx.data[x][y]), 3):
          print("{:e} {:e} {:e} ".format(dx.data[x][y][z], dx.data[x][y][z+1], dx.data[x][y][z+2]), file=f)
    print(*dx.file_footer, sep='', end='', file=f)
  finally:
    if filename:
      f.close()

def readDX(filename):
  data = []
  origin = None
  delta  = None
  file_header=[]
  file_footer=[]

  with open(filename) as f:
    isData = False
    datapoints = 0
    x, y, z = (0, 0, 0)
    for line in f:
      if (not isData):
        file_header.append(line)
      elif (isData and datapoints <= 0):
        file_footer.append(line)

      if isData and datapoints > 0:
        for val in line.split():
          data[x][y][z] = float(val)
          datapoints -= 1

          if datapoints == 0:
            continue

          # increase coordinates
          z += 1
          if z >= len(data[x][y]):
            z = 0
            y += 1

          if y >= len(data[x]):
            y = 0
            x += 1

          if x >= len(data):
            raise Exception('x too large')
      elif not isData:
        if line.startswith('origin ') or line.startswith('delta '):
          m = re.search('^(?P<type>origin|delta)\s+(?P<x>\S+)\s+(?P<y>\S+)\s+(?P<z>\S+)$', line)
          if m:
            p = Point(float(m.group('x')), float(m.group('y')), float(m.group('z')))

            if m.group('type') == 'origin':
              origin = p
            elif m.group('type') == 'delta':
              # if delta exist, add new values to existing
              if delta:
                delta = delta.add(p)
              else:
                delta = p

        elif line.startswith('object 2'):
          # array dimensions
          p = re.compile(' gridconnections counts (?P<x>\d+) (?P<y>\d+) (?P<z>\d+)$')

          m = p.search(line)
          if m:
            data = [ [ [7 for i in range(int(m.group('z'))) ] for j in range(int(m.group('y'))) ] for k in range(int(m.group('x'))) ]
            datapoints = 1
            for coord in ('x', 'y', 'z'):
              datapoints *= int(m.group(coord))
        elif line.rstrip().endswith(' data follows'):
          isData = True

  return OpenDX(data = data, origin = origin, delta = delta, file_header = file_header, file_footer = file_footer)

def electrostatic_hull(data, iterations = 1, filename = ''):
  newdata = data
  for i in range(iterations):
    data = copy.deepcopy(newdata)
    res = set()
    for x in range(len(data)):
      for y in range(len(data[x])):
        for z in range(len(data[x][y])):
          if data[x][y][z] == 0:
            # look at datapoints around this one

            # x-axis
            for nx in (x-1,x+1):
              if nx >= 0 and nx < len(data) and data[nx][y][z] > 0:
                newdata[nx][y][z] = 0
                res.add((nx, y, z))

            # y-axis
            for ny in (y-1,y+1):
              if ny >= 0 and ny < len(data) and data[x][ny][z] > 0:
                res.add((x, ny, z))
                newdata[x][ny][z] = 0

            # z-axis
            for nz in (z-1,z+1):
              if nz >= 0 and nz < len(data) and data[x][y][nz] > 0:
                res.add((x, y, nz))
                newdata[x][y][nz] = 0

    if i == iterations - 1:
      f = sys.stdout
      if filename:
        f = open(filename, 'w')

      try:
        # print csv header
        print("x,y,z", file=f)
        for t in sorted(res):
          print("{},{},{}".format(*t), file=f)
      finally:
        if filename:
          f.close()

def dx2csv(dx_list, filter = None, ids = None, sep = ' ', filename = ''):
  f = sys.stdout
  if filename:
    f = open(filename, 'w')

  if filter:
    # make sure the points are sorted and unique
    filter = sorted(set(filter))

  if type(dx_list) == OpenDX:
    dx_list = [dx_list]

  try:
    # print header
    if len(dx_list) > 0:
      dx = dx_list[0]

      # print empty column for row names
      # print(sep, end='', file=f)

      if filter:
        for p in filter:
          print("{}{}.{}.{}".format(sep, p.x, p.y, p.z), end='', file=f)
      else:
        for x in range(len(dx.data)):
          for y in range(len(dx.data[x])):
            for z in range(len(dx.data[x][y])):
              print("{}{}.{}.{}".format(sep, x, y, z), end='', file=f)

      # print newline
      print(file=f)

    # print data
    for i in range(len(dx_list)):
      dx = dx_list[i]
      if ids:
        print(ids[i], end='', file=f)

      if filter:
        for p in filter:
          print("{}{:e}".format(sep, dx.data[p.x][p.y][p.z]), end='', file=f)
      else:
        for x in range(len(dx.data)):
          for y in range(len(dx.data[x])):
            for z in range(len(dx.data[x][y])):
              print("{}{:e}".format(sep, dx.data[x][y][z]), end='', file=f)

      # print newline
      print(file=f)

  finally:
    if filename:
      f.close()

def csv2points(filename):
  points = []
  with open(filename) as f:
    for line in f:
      if points or re.match('[0-9]', line):
        points.append(Point(*map(int, line.rstrip().split(','))))

  return points

# main
# try:
#   snakemake
# except NameError:
#   snakemake = None
#
# if snakemake:
#   if hasattr(snakemake.params, 'mode'):
#     mode = snakemake.params['mode']
#     if mode == 'calculate_points':
#       distance = snakemake.wildcards['distance']
#       dx = readDX(snakemake.input['dx_smol'])
#       electrostatic_hull(data = dx.data, iterations = int(int(distance) / dx.delta[0]) + 1 , filename = snakemake.output[0])
#     elif mode == 'generate_data_csv':
#       filter = None
#       if hasattr(snakemake.input, 'eh_points'):
#         filter = csv2points(snakemake.input['eh_points'])
#
#       ids = []
#       dx_list = []
#       p = re.compile('/(?:pot|smol)_(\S+)\.dx$')
#       for f in [snakemake.input['dx']] if type(snakemake.input['dx']) == str else snakemake.input['dx']:
#         dx_list.append(readDX(f))
#         ids.append(p.search(f).group(1))
#
#       dx2csv(dx_list, filter = filter, ids = ids, filename = snakemake.output['csv'])
