#' Rotating calipers algorithm
#'
#' @description
#' credits go to Daniel Wollschlaeger <https://github.com/ramnathv>
#'
#' @importFrom grDevices chull
getMinBBox <- function(xy) {
  stopifnot(is.matrix(xy), is.numeric(xy), nrow(xy) >= 2, ncol(xy) == 2)

  ## rotating calipers algorithm using the convex hull
  H    <- grDevices::chull(xy)      ## hull indices, vertices ordered clockwise
  n    <- length(H)      ## number of hull vertices
  hull <- xy[H, ]        ## hull vertices

  ## unit basis vectors for all subspaces spanned by the hull edges
  hDir  <- diff(rbind(hull, hull[1, ])) ## hull vertices are circular
  hLens <- sqrt(rowSums(hDir^2))        ## length of basis vectors
  huDir <- diag(1/hLens) %*% hDir       ## scaled to unit length

  ## unit basis vectors for the orthogonal subspaces
  ## rotation by 90 deg -> y' = x, x' = -y
  ouDir <- cbind(-huDir[ , 2], huDir[ , 1])

  ## project hull vertices on the subspaces spanned by the hull edges, and on
  ## the subspaces spanned by their orthogonal complements - in subspace coords
  projMat <- rbind(huDir, ouDir) %*% t(hull)

  ## range of projections and corresponding width/height of bounding rectangle
  rangeH  <- matrix(numeric(n*2), ncol=2)  ## hull edge
  rangeO  <- matrix(numeric(n*2), ncol=2)  ## orthogonal subspace
  widths  <- numeric(n)
  heights <- numeric(n)

  for(i in seq(along=numeric(n))) {
    rangeH[i, ] <- range(projMat[  i, ])

    ## the orthogonal subspace is in the 2nd half of the matrix
    rangeO[i, ] <- range(projMat[n+i, ])
    widths[i]   <- abs(diff(rangeH[i, ]))
    heights[i]  <- abs(diff(rangeO[i, ]))
  }

  ## extreme projections for min-area rect in subspace coordinates
  ## hull edge leading to minimum-area
  eMin  <- which.min(widths*heights)
  hProj <- rbind(   rangeH[eMin, ], 0)
  oProj <- rbind(0, rangeO[eMin, ])

  ## move projections to rectangle corners
  hPts <- sweep(hProj, 1, oProj[ , 1], "+")
  oPts <- sweep(hProj, 1, oProj[ , 2], "+")

  ## corners in standard coordinates, rows = x,y, columns = corners
  ## in combined (4x2)-matrix: reverse point order to be usable in polygon()
  ## basis formed by hull edge and orthogonal subspace
  basis <- cbind(huDir[eMin, ], ouDir[eMin, ])
  hCorn <- basis %*% hPts
  oCorn <- basis %*% oPts
  pts   <- t(cbind(hCorn, oCorn[ , c(2, 1)]))

  ## angle of longer edge pointing up
  dPts <- diff(pts)
  e    <- dPts[which.max(rowSums(dPts^2)), ] ## one of the longer edges
  eUp  <- e * sign(e[2])       ## rotate upwards 180 deg if necessary
  deg  <- atan2(eUp[2], eUp[1])*180 / pi     ## angle in degrees

  return(list(pts=pts, width=heights[eMin], height=widths[eMin], angle=deg))
}


#' @importFrom grDevices chull
getBBoxAngle = function(vertices, alpha) {
  centroid = apply(vertices, 2, mean)

  angleRadians = (90-alpha) * pi/180

  rotatedX = cos(angleRadians) * (vertices[, 1] - centroid[1]) -
             sin(angleRadians) * (vertices[, 2] - centroid[2]) +
             centroid[1]
  rotatedY = sin(angleRadians) * (vertices[, 1] - centroid[1]) +
             cos(angleRadians) * (vertices[, 2] - centroid[2]) +
             centroid[2]

  bBox = cbind(range(rotatedX), range(rotatedY))
  height = diff(bBox[,2])
  width = diff(bBox[,1])

  bBox = rbind(bBox[1,],
               c(bBox[1, 1], bBox[2, 2]),
               (bBox[2, ]),
               c(bBox[2, 1], bBox[1, 2]))

  bBoxUnrotatedX = cos(-angleRadians) * (bBox[, 1] - centroid[1]) -
                   sin(-angleRadians) * (bBox[, 2] - centroid[2]) +
                   centroid[1]
  bBoxUnrotatedY = sin(-angleRadians) * (bBox[, 1] - centroid[1]) +
                   cos(-angleRadians) * (bBox[, 2] - centroid[2]) +
                   centroid[2]
  unrotatedBbox = cbind(bBoxUnrotatedX, bBoxUnrotatedY)
  # plot(unrotatedBbox, pch=3)
  # polygon(unrotatedBbox)
  # points(vertices)
  # polygon(minBbox$pts, border='yellow', lty='dashed')
  list(angle=alpha, height=height, width=width, pts=unrotatedBbox)
}


adjustAcuteAngles = function(xy, angle, minAngle) {
  xy_mat = as.matrix(xy[,1:2])
  rads = angle*pi/180
  nPoints = nrow(xy)

  # Using cosine rule on vectors
  # acos(theta) = (AB_vec x BC_vec) / (|AB| * |BC|)
  vectors = xy[-1,]-xy[-nPoints,]

  # Final angles for each point
  angles = getAngles(xy_mat)

  # Mask angles less than minimum
  mask = angles <= minAngle
  mask[is.na(mask)] = FALSE
  # Index of point with angles less than minimum
  indices = (seq_along(mask)+1)[mask]

  indicesOffset = (indices %% 2 == 1) * 2 - 1

  # Calculate perpendicular lines for points
  x = xy[indices, 1]
  y = xy[indices, 2]
  nIndices = length(indices)
  a = rep(tan(rads), nIndices)
  b = y - a * x
  a_perpendicular = -1/a
  b_perpendicular = y - a_perpendicular*x


  b2 = b_perpendicular
  a2 = a_perpendicular

  # Intersect previous straight line with the perpendicular
  x = xy[indices-indicesOffset, 1]
  y = xy[indices-indicesOffset, 2]
  a1 = a
  b1 = y - a1 * x
  xs = (b2-b1)/(a1-a2)
  ys = a2*xs+b2

  # Replace angled points
  xy[indices-indicesOffset, 1] = xs
  xy[indices-indicesOffset, 2] = ys
  return(xy)
}

outerCurvePoints = function(waypoints, angle, flightLineDistance) {
  mask = getAngles(waypoints) == 180
  mask[is.na(mask)] = TRUE

  rads = angle*pi/180
  angularCoef = tan(rads)
  nWaypoints = nrow(waypoints)

  oddPoints = seq(1, nWaypoints, 2)
  evenPoints = seq(2, nWaypoints, 2)

  offsetX = waypoints[evenPoints,1]-waypoints[oddPoints,1]

  intercept = waypoints[oddPoints, 2] - angularCoef*waypoints[oddPoints, 1]
  xMove = rep(flightLineDistance / sqrt(1 + angularCoef**2), nWaypoints)
  isNegative = rep(offsetX < 0, each=2)
  isNegative[seq(1, nWaypoints, 2)] = !isNegative[seq(1, nWaypoints, 2)]
  xMove[isNegative] = -xMove[isNegative]
  yMove = xMove * angularCoef

  xs = waypoints[, 1] + xMove
  ys = waypoints[, 2] + yMove

  curvePoints = as.data.frame(cbind(xs, ys))
  curvePoints$index = seq(1, nWaypoints)
  curvePoints$before = (curvePoints$index %% 2) == 1
  curvePoints = curvePoints[!c(FALSE, mask),]

  return(curvePoints)
}

getAngles = function(waypoints) {
  nPoints = nrow(waypoints)
  vectors = waypoints[-1,]-waypoints[-nPoints,]
  dists = sqrt(apply(vectors ** 2, 1, sum))

  ab = vectors[-nPoints+1,]
  bc = vectors[-1,]
  dist_ab = dists[-nPoints+1]
  dist_bc = dists[-1]
  dotProds = sapply(seq_len(nPoints-2), function(x) ab[x,] %*% bc[x,])
  # Final angles for each point
  angles = suppressWarnings(180-round(acos(dotProds/(dist_ab*dist_bc))*180/pi,2))
  angles
}
