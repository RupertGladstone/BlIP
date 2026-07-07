
yearinsec = 365.25*24.0*60.0*60.0
rhoi = 910.0/(1.0e6*yearinsec^2)
rhow = 1027.0/(1.0e6*yearinsec^2)

gravity = -9.81*yearinsec^2

GLTolerance = 1.0e-4

mm = 1.0/3.0
nn = 3.0

MINH = 10.0

beta = 0.1/1000.0

eta = (1.0*100.0)^(-1.0/nn)
--! approx minus 5C

-- MISMIP+ bedrock elevation, Eq. (1) and Table 1 of Cornford et al. (2020)
xbar = 300.0e3
B0   = -150.0
B2   = -728.8
B4   = 343.91
B6   = -50.75
wc   = 24.0e3
fc   = 4.0e3
dc   = 500.0

function mismip_bedrock(x, y)
  local Bx = B0 + B2*(x/xbar)^2 + B4*(x/xbar)^4 + B6*(x/xbar)^6
  local By = dc*( 1.0/(1.0 + math.exp(-2.0*(y - wc)/fc))
                + 1.0/(1.0 + math.exp( 2.0*(y + wc)/fc)) )
  return math.max(Bx + By, -720.0)
end

function Zb(x, y, H)
  local base_afloat = -(rhoi/rhow) * H
  return math.max(mismip_bedrock(x, y), base_afloat)
end

function Zs(x, y, H)
  return Zb(x, y, H) + H
end
