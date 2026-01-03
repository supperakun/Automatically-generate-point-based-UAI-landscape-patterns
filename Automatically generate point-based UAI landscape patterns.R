###Custom functions

####five functions used to caculate tree-based UAI
Near.f=function(x0,y0,data1){
  x1=abs(data1$x-x0); y1=abs(data1$y-y0)
  dis=sqrt(x1*x1+y1*y1)
  nord=order(dis)
  nord[1:5]
  data2=rbind(data1[nord[1],],data1[nord[2],],data1[nord[3],],data1[nord[4],],data1[nord[5],])
}
Angle.f=function(x0,y0,xi,yi){
  deltax=xi-x0
  deltay=yi-y0
  if (deltay>0 &deltax>=0) {Angle=atan(deltax/deltay)} else
    if (deltay==0 &deltax>=0) {Angle=pi*0.5} else
      if (deltay<=0 &deltax>=0) {Angle=pi+atan(deltax/deltay)} else
        if (deltay<0 &deltax<0) {Angle=pi+atan(deltax/deltay)} else
          if (deltay==0 &deltax<0) {Angle=pi*1.5} else
          {Angle=pi*2+atan(deltax/deltay)}
  Angle
}
AngleDiff.f=function(Angle1,Angle2){
  AngleDiff=Angle1-Angle2
  if (AngleDiff>pi) {AngleDiff=2*pi-AngleDiff}
  if (AngleDiff>=72*pi/180){AngleDiff=0} else
  {AngleDiff=1}
}
w.f=function(x0,y0,x,y){
  Angle4=c()
  for(i in 1:4){
    xi=x[i];yi=y[i]
    Angle4[i]=Angle.f(x0,y0,xi,yi)
  }
  Angle4o=order(Angle4)
  count12=AngleDiff.f(Angle4[Angle4o[2]],Angle4[Angle4o[1]])
  count23=AngleDiff.f(Angle4[Angle4o[3]],Angle4[Angle4o[2]])
  count34=AngleDiff.f(Angle4[Angle4o[4]],Angle4[Angle4o[3]])
  count41=AngleDiff.f(Angle4[Angle4o[4]],Angle4[Angle4o[1]])
  count=(count12+count23+count34+count41)/4
  w=count
}
W.even=function(x0,y0,Near4){
  for(i in 1:4){
    xi=Near4$x[i];yi=Near4$y[i]
    Near4$angle[i]=Angle.f(x0,y0,xi,yi)
  }
  Near4o=Near4[order(Near4[,"angle"]),]
  Near4o$count=NA
  for (i in 1:3){
    Near4o$count[i]=AngleDiff.f(Near4o$angle[i+1],Near4o$angle[i])
  }
  Near4o$count[4]=AngleDiff.f(Near4o$angle[4],Near4o$angle[1])
  evenNo=c(Near4o[Near4o$count==0,]$No)
  return(evenNo)
}
W.cluster=function(x0,y0,Near4){
  for(i in 1:4){
    xi=Near4$x[i];yi=Near4$y[i]
    Near4$angle[i]=Angle.f(x0,y0,xi,yi)
  }
  Near4o=Near4[order(Near4[,"angle"]),]
  Near4o$count=NA
  for (i in 1:3){
    Near4o$count[i]=AngleDiff.f(Near4o$angle[i+1],Near4o$angle[i])
  }
  Near4o$count[4]=AngleDiff.f(Near4o$angle[4],Near4o$angle[1])
  clusterNo=c(Near4o[Near4o$count==1,]$No)
  return(clusterNo)
}
##########end of tree-based UAI function

####functions used to caculate double-arc
circle=function(x,y){
  point1=x
  point2=y
  known.pair1 <- structure(c(point1$x, point2$x, 
                             point1$y, point2$y), .Dim = c(2L, 2L), 
                           .Dimnames = list(NULL, c("x", "y")))
  ymin <- min(known.pair1[,2])
  ymax <- max(known.pair1[,2])
  xmin <- min(known.pair1[,1])
  xmax <- max(known.pair1[,1])
  if(known.pair1[1,2]==known.pair1[2,2]){
    known.pair=known.pair1
  }else{
    known.pair[1,]=known.pair1[known.pair1[,2]==ymin,]
    known.pair[2,]=known.pair1[known.pair1[,2]==ymax,]}
  dif = diff(known.pair)
  ## Distance and angle (/_KkB) between the two known points
  d1 <- sqrt(sum(dif^2))
  r <- d1/(2*sin(0.4*pi))
  theta2 <- pi/10
  if(dif[2]==0){
    ## Find center of one circle (using /_BkC1)
    dx1 <- cos(theta2)*r
    dy1 <- sin(theta2)*r
    p1 <- known.pair[known.pair[,1]==xmin,] + c(dx1, dy1)
    p1 <- st_point(c(p1[1],p1[2]),dim = "XY")
    c1 <- st_buffer(p1, r,3000)
    ## Find center of other circle (using /_BkC2)
    dx2 <- dx1
    dy2 <- -dy1
    p2 <- known.pair[known.pair[,1]==xmin,] + c(dx2, dy2)
    p2 <- st_point(c(p2[1],p2[2]),dim = "XY")
    c2 <- st_buffer(p2, r,3000)
    c <- st_union(c1, c2)
    c
  }else if(dif[1]>0){
    theta1 <- atan(do.call("/", as.list(rev(abs(dif)))))
    ## Find center of one circle (using /_BkC1)
    dx1 <- cos(theta1+theta2)*r
    dy1 <- sin(theta1+theta2)*r
    p1 <- known.pair[1,] + c(dx1, dy1)
    p1 <- st_point(c(p1[1],p1[2]),dim = "XY")
    c1 <- st_buffer(p1, r,3000)
    ## Find center of other circle (using /_BkC2)
    dx2 <- cos(theta1-theta2)*r
    dy2 <- sin(theta1-theta2)*r
    p2 <- known.pair[1,] + c(dx2, dy2)
    p2 <- st_point(c(p2[1],p2[2]),dim = "XY")
    c2 <- st_buffer(p2, r,3000)
    c <- st_union(c1, c2)
    c
  }else if(dif[1]<0){
    theta1 <- atan(do.call("/", as.list(rev(abs(dif)))))
    ## Find center of one circle (using /_BkC1)
    dx1 <- -cos(theta1+theta2)*r
    dy1 <- sin(theta1+theta2)*r
    p1 <- known.pair[1,] + c(dx1, dy1)
    p1 <- st_point(c(p1[1],p1[2]),dim = "XY")
    c1 <- st_buffer(p1, r,3000)
    ## Find center of other circle (using /_BkC2)
    dx2 <- -cos(theta1-theta2)*r
    dy2 <- sin(theta1-theta2)*r
    p2 <- known.pair[1,] + c(dx2, dy2)
    p2 <- st_point(c(p2[1],p2[2]),dim = "XY")
    c2 <- st_buffer(p2, r,3000)
    c <- st_union(c1, c2)
    c
  }else if(dif[1]==0){
    ## Find center of one circle (using /_BkC1)
    dx1 <- sin(theta2)*r
    dy1 <- cos(theta2)*r
    p1 <- known.pair[1,] + c(dx1, dy1)
    p1 <- st_point(c(p1[1],p1[2]),dim = "XY")
    c1 <- st_buffer(p1, r,3000)
    ## Find center of other circle (using /_BkC2)
    dx2 <- -dx1
    dy2 <- dy1
    p2 <- known.pair[1,] + c(dx2, dy2)
    p2 <- st_point(c(p2[1],p2[2]),dim = "XY")
    c2 <- st_buffer(p2, r,3000)
    c <- st_union(c1, c2)
    c
  }
}

#######functions used to point-based UAI 
st_ends_heading <- function(line){
  M <- sf::st_coordinates(line)
  i <- c(2, nrow(M) - 1)
  j <- c(1, -1)
  
  headings <- mapply(i, j, FUN = function(i, j) {
    Ax <- M[i-j,1]
    Ay <- M[i-j,2]
    Bx <- M[i,1]
    By <- M[i,2]
    unname(atan2(Ay-By, Ax-Bx))
  })
  
  return(headings)
}

st_extend_line <- function(line, distance, end = "BOTH"){
  if (!(end %in% c("BOTH", "HEAD", "TAIL")) | length(end) != 1) stop("'end' must be 'BOTH', 'HEAD' or 'TAIL'")
  
  M <- sf::st_coordinates(line)[,1:2]
  keep <- end == c("TAIL", "HEAD")
  ends <- c(1, nrow(M))[keep]
  headings <- st_ends_heading(line)[keep]
  distances <- if (length(distance) == 1) rep(distance, 2) else rev(distance[1:2])
  
  M[ends,] <- M[ends,] + distances[keep] * c(cos(headings), sin(headings))
  newline <- sf::st_linestring(M)
  
  # If input is sfc_LINESTRING and not sfg_LINESTRING
  if (is.list(line)) newline <- sf::st_sfc(newline, crs = sf::st_crs(line))
  
  return(newline)
}

ptincerts2V=function(pintpair){
  deltax = pintpair[1,1]-pintpair[2,1]
  deltay = pintpair[1,2]-pintpair[2,2]
  k=deltay/deltax
  b=pintpair[1,2]-(pintpair[1,1]*k)
  if(abs(deltax)<0.00001){
    LP=c(pintpair[1,1],-1000)
    RP=c(pintpair[1,1],1000)
  }
  if(abs(deltay)<0.00001){
    LP=c(-1000,pintpair[1,2])
    RP=c(1000,pintpair[1,2])
  }
  if(abs(deltay)>0.0001&abs(deltax)>0.0001){
    if(k>0){
      RP=c(1000,1000*k+b)
      LP=c(-1000,-1000*k+b)
    }else{
      RP=c(-1000,-1000*k+b)
      LP=c(1000,1000*k+b)
    }
  }
  ptincerts=rbind(LP,RP)
}

ptincerts2H=function(pintpair){
  deltax = pintpair[1,1]-pintpair[2,1]
  deltay = pintpair[1,2]-pintpair[2,2]
  k=deltay/deltax
  b=pintpair[1,2]-(pintpair[1,1]*k)
  if(abs(deltax)<0.00001){
    LP=c(pintpair[1,1],-1000)
    RP=c(pintpair[1,1],1000)
  }
  if(abs(deltay)<0.00001){
    LP=c(-1000,pintpair[1,2])
    RP=c(1000,pintpair[1,2])
  }
  if(abs(deltay)>0.0001&abs(deltax)>0.0001){
    if(pintpair[2,2]<pintpair[1,2]){
      if(k>0){
        RP=c(1000,1000*k+b)
        LP=c(-1000,-1000*k+b)
      }else{
        RP=c(-1000,-1000*k+b)
        LP=c(1000,1000*k+b)
      }
    }else if(pintpair[2,2]>pintpair[1,2]){
      if(k<0){
        RP=c(1000,1000*k+b)
        LP=c(-1000,-1000*k+b)
      }else{
        RP=c(-1000,-1000*k+b)
        LP=c(1000,1000*k+b)
      }
    }
  }
  ptincerts=rbind(LP,RP)
}

ptincerts1L=function(pintpair){
  deltax = pintpair[1,1]-pintpair[2,1]
  deltay = pintpair[1,2]-pintpair[2,2]
  k=deltay/deltax
  b=pintpair[1,2]-(pintpair[1,1]*k)
  if(abs(deltax)<0.00001){
    LP=c(pintpair[1,1],1000)
    RP=c(pintpair[1,1],-1000)
  }
  if(abs(deltay)<0.00001){
    LP=c(-1000,pintpair[1,2])
    RP=c(1000,pintpair[1,2])
  }
  if(abs(deltay)>0.0001&abs(deltax)>0.0001){
    if(pintpair[2,2]<pintpair[1,2]){
      if(k>0){
        RP=c(1000,1000*k+b)
        LP=c(-1000,-1000*k+b)
      }else{
        RP=c(-1000,-1000*k+b)
        LP=c(1000,1000*k+b)
      }
    }else if(pintpair[2,2]>pintpair[1,2]){
      if(k<0){
        RP=c(1000,1000*k+b)
        LP=c(-1000,-1000*k+b)
      }else{
        RP=c(-1000,-1000*k+b)
        LP=c(1000,1000*k+b)
      }
    }
  }
  ptincerts=rbind(LP,RP)
}

ptincerts1R=function(pintpair){
  deltax = pintpair[1,1]-pintpair[2,1]
  deltay = pintpair[1,2]-pintpair[2,2]
  k=deltay/deltax
  b=pintpair[1,2]-(pintpair[1,1]*k)
  if(abs(deltax)<0.00001){
    LP=c(pintpair[1,1],1000)
    RP=c(pintpair[1,1],-1000)
  }
  if(abs(deltay)<0.00001){
    LP=c(-1000,pintpair[1,2])
    RP=c(1000,pintpair[1,2])
  }
  if(abs(deltay)>0.0001&abs(deltax)>0.0001){
    RP=c(1000,1000*k+b)
    LP=c(-1000,-1000*k+b)
  }
  ptincerts=rbind(LP,RP)
}



ptincerts3V=function(pintpair1,pintpair2){
  fit1 <- lm(y~x,data = pintpair1)
  fit2 <- lm(y~x,data = pintpair2)
  coef1 = coef(fit1)
  coef2 = coef(fit2)
  coef1[is.na(coef1)] <- 0
  coef2[is.na(coef2)] <- 0
  if(abs(coef1[2]-coef2[2])<0.0000001){
    if(abs(coef1[2])<0.0000001){
      RP1=c(pintpair1[2,1],-1000)
      LP1=c(pintpair1[1,1],1000)
      RP2=c(pintpair1[2,1],1000)
      LP2=c(pintpair1[1,1],-1000)
      ptincerts=rbind(RP1,LP1,RP2,LP2)}
    else{
      RP1=c(1000,1000*coef1[2]+coef1[1])
      LP1=c(-1000,-1000*coef1[2]+coef1[1])
      RP2=c(1000,1000*coef2[2]+coef2[1])
      LP2=c(-1000,-1000*coef2[2]+coef2[1])
      ptincerts=rbind(RP1,LP1,RP2,LP2)}  
  }else{
    A <- matrix(c(coef1[2], -1,
                  coef2[2], -1), byrow = T, nrow = 2)
    b <- c(-coef1[1], -coef2[1])
    MP =  solve(A, b)
    if(coef1[2]>0){
      LP1=c(1000,1000*coef1[2]+coef1[1])
      RP1=c(-1000,-1000*coef1[2]+coef1[1])
    }
    if(coef1[2]==0){
      LP1=c(1000,pintpair1[1,2])
      RP1=c(-1000,pintpair1[1,2])
      if(is.na(coef(fit1)[2])){
        MP =c(pintpair1[1,1],pintpair1[1,1]*coef2[2]+coef2[1])
        LP1=c(pintpair1[1,1],1000)
        RP1=c(pintpair1[1,1],-1000)
      }
    }
    if(coef1[2]<0){
      RP1=c(1000,1000*coef1[2]+coef1[1])
      LP1=c(-1000,-1000*coef1[2]+coef1[1])
    }
    if(coef2[2]>0){
      LP2=c(1000,1000*coef2[2]+coef2[1])
      RP2=c(-1000,-1000*coef2[2]+coef2[1])
    }
    if(coef2[2]==0){
      LP2=c(1000,pintpair2[1,2])
      RP2=c(-1000,pintpair2[1,2])
      if(is.na(coef(fit2)[2])){
        MP =c(pintpair2[1,1],pintpair2[1,1]*coef1[2]+coef1[1])
        LP2=c(pintpair2[1,1],1000)
        RP2=c(pintpair2[1,1],-1000)
      }
    }
    if(coef2[2]<0){
      RP2=c(1000,1000*coef2[2]+coef2[1])
      LP2=c(-1000,-1000*coef2[2]+coef2[1])
    }
    ptincerts=rbind(MP,RP1,LP1,RP2,LP2)
  }
}

ptincerts3H=function(pintpair1,pintpair2){
  fit1 <- lm(y~x,data = pintpair1)
  fit2 <- lm(y~x,data = pintpair2)
  coef1 = coef(fit1)
  coef2 = coef(fit2)
  coef1[is.na(coef1)] <- 0
  coef2[is.na(coef2)] <- 0
  if(abs(coef1[2]-coef2[2])<0.0000001){
    RP1=c(1000,1000*coef1[2]+coef1[1])
    LP1=c(-1000,-1000*coef1[2]+coef1[1])
    RP2=c(1000,1000*coef2[2]+coef2[1])
    LP2=c(-1000,-1000*coef2[2]+coef2[1])
    ptincerts=rbind(RP1,LP1,RP2,LP2)  
  }else{
    A <- matrix(c(coef1[2], -1,
                  coef2[2], -1), byrow = T, nrow = 2)
    b <- c(-coef1[1], -coef2[1])
    MP =  solve(A, b)
    RP1=c(1000,1000*coef1[2]+coef1[1])
    LP1=c(-1000,-1000*coef1[2]+coef1[1])
    RP2=c(1000,1000*coef2[2]+coef2[1])
    LP2=c(-1000,-1000*coef2[2]+coef2[1]) 
    ptincerts=rbind(MP,RP1,LP1,RP2,LP2)
  }
}

typecheck=function(pintpair1,pintpair2){
  fit1 <- lm(y~x,data = pintpair1)
  fit2 <- lm(y~x,data = pintpair2)
  coef1 = coef(fit1)
  coef2 = coef(fit2)
  coef1[is.na(coef1)] <- 0
  coef2[is.na(coef2)] <- 0
  RUCHECK=pintpair2[1,1]*coef1[2]+coef1[1]
  LDCHECK=pintpair2[2,1]*coef1[2]+coef1[1]
  LUCHECK=(pintpair1[1,2]-coef2[1])/coef2[2]
  RDCHECK=(pintpair1[2,2]-coef2[1])/coef2[2]
  if(coef1[2]>0){
    RUJUG=RUCHECK<pintpair2[1,2]
    LDJUG=LDCHECK>pintpair2[2,2]
  }else if(coef1[2]<0){
    RUJUG=RUCHECK<pintpair2[1,2]
    LDJUG=LDCHECK>pintpair2[2,2]
  }
  if(abs(coef1[2])<0.00000001){
    RUCHECK=pintpair1[1,1]
    LDCHECK=pintpair1[2,1]
    if(is.na(coef(fit1)[2])){
      RUJUG=pintpair2[1,1]>pintpair1[1,1]
      LDJUG=pintpair2[2,1]<pintpair1[2,1] 
    }else{
      RUJUG=RUCHECK<pintpair2[1,1]
      LDJUG=LDCHECK>pintpair2[2,1]}
  }
  if(coef2[2]>0){
    LUJUG=LUCHECK>pintpair1[1,1]
    RDJUG=RDCHECK<pintpair1[2,1]
  }else if(coef2[2]<0){
    LUJUG=LUCHECK>pintpair1[1,1]
    RDJUG=RDCHECK<pintpair1[2,1]
  }
  if(abs(coef2[2])<0.00000001){
    LUCHECK=pintpair2[1,1]
    RDCHECK=pintpair2[2,1]
    if(is.na(coef(fit2)[2])){
      LUJUG=pintpair1[1,1]<pintpair2[1,1]
      RDJUG=pintpair1[2,1]>pintpair2[2,1] 
    }else{
      LUJUG=pintpair1[1,2]>pintpair1[1,2]
      RDJUG=pintpair1[2,2]<pintpair1[2,2]}
  }
  JUG=c(RUJUG,LDJUG,LUJUG,RDJUG)
}

Anglexaxis=function(xi,yi){
  if (yi>0 &xi>=0) {Angle=atan(yi/xi)/pi*180} else
    if (abs(yi)<0.0000001 &xi>=0) {Angle=0} else
      if (yi<=0 &xi>=0) {Angle=360+atan(yi/xi)/pi*180} else
        if (yi<0 &xi<0) {Angle=180+atan(yi/xi)/pi*180} else
          if (abs(yi)<0.0000001 &xi<0) {Angle=180} else 
            if (yi<0 &abs(xi)<0.0000001) {Angle=270} else
              if (yi>0 &abs(xi)<0.0000001) {Angle=90} else
              {Angle=180+atan(yi/xi)/pi*180}
  Angle
}

ptV=function(pintpair){
  deltax = pintpair[1,1]-pintpair[2,1]
  deltay = pintpair[1,2]-pintpair[2,2]
  k=deltay/deltax
  b=pintpair[1,2]-(pintpair[1,1]*k)
  if(abs(deltax)<0.00001){
    UP=c(pintpair[1,1],1000)
    DP=c(pintpair[1,1],-1000)
  }else{
    if(abs(deltax)>=abs(deltay)){
      k=deltay/deltax
      b=pintpair[1,2]-(pintpair[1,1]*k)
      LP=c(-1000,-1000*k+b)
      RP=c(1000,1000*k+b)
      ptincerts=rbind(LP,RP)
    }else{
      if(k>0){
        UP=c(1000,1000*k+b)
        DP=c(-1000,-1000*k+b) 
      }else{
        DP=c(1000,1000*k+b)
        UP=c(-1000,-1000*k+b)
      }
      ptincerts=rbind(UP,DP) 
    }
    
  }
  
}

extend_point <- function(A, B, distance = 10000) {
  
  # 计算向量AB
  AB <- B - A
  
  # 计算向量AB的模（长度）
  norm_AB <- sqrt(sum(AB^2))
  
  # 检查是否A和B重合
  if (norm_AB == 0) {
    stop("Points A and B are identical, cannot determine direction.")
  }
  
  # 计算反向单位向量
  unit_vector <- -AB / norm_AB
  
  # 计算点C的坐标：A + distance * (-单位向量)
  C <- A + distance * unit_vector
  
  return(C)
}

is_convex_quadrilateral <- function(points) {
  extended_points <- rbind(points, points[1, ],points[2, ])
  
  # 初始化叉积符号向量
  cross_signs <- numeric(4)
  
  # 遍历每个顶点计算转向
  for (i in 1:4) {
    # 获取连续的三个点
    a <- extended_points[i, ]
    b <- extended_points[i+1, ]
    c <- extended_points[i+2, ]
    
    # 计算向量AB和BC
    ab <- c(b$x - a$x, b$y - a$y)
    bc <- c(c$x - b$x, c$y - b$y)
    
    # 计算二维叉积（z分量）
    cross_product <- ab[1] * bc[2] - ab[2] * bc[1]
    
    # 记录符号（保留原始数值用于调试）
    cross_signs[i] <- sign(cross_product)
  }
  
  # 判断所有符号是否一致
  if (all(cross_signs >= 0)) {
    return(TRUE)
  } else if (all(cross_signs <= 0)) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}



move_points_to_centroid <- function(points) {
  
  centroid <- c(median(points[,1]),median(points[,2]))
  
  # 计算每个点到重心的欧几里得距离
  distances <- sqrt(rowSums((points - centroid)^2))
  
  # 找到最小距离
  min_dist <- min(distances)
  
  # 初始化新点坐标矩阵
  new_points <- points
  
  # 对每个点进行移动
  for (i in 1:4) {
    if (distances[i] > 0) {  # 仅对不在重心的点进行移动
      # 计算单位方向向量：从点指向重心
      direction <- (centroid - points[i, ]) / distances[i]
      # 新位置 = 原位置 + 最小距离 * 方向向量
      new_points[i, ] <- points[i, ] + min_dist * direction
    }
  }
  
  # 返回新点的坐标
  return(new_points)
}

sort_points_clockwise <- function(df) {
   
  # 计算几何中心
  centroid <- c(mean(df$x), mean(df$y))
  
  # 计算相对极角（弧度）
  calc_angle <- function(point) {
    dx <- as.numeric(point[1]) - centroid[1]
    dy <- as.numeric(point[2]) - centroid[2]
    atan2(dy, dx)  # 计算结果范围：-π到π
  }
  
  # 添加临时极角列
  df$polar_angle <- apply(df, 1, calc_angle)
  
  # 按极角逆时针排序（默认排序方向）
  ordered_df <- df[order(-df$polar_angle), ]  # 负号实现降序排列
  
  # 计算有符号面积验证方向
  signed_area <- 0
  n <- nrow(ordered_df)
  for (i in 1:n) {
    x_i <- ordered_df$x[i]
    y_i <- ordered_df$y[i]
    x_next <- ordered_df$x[i %% n + 1]
    y_next <- ordered_df$y[i %% n + 1]
    signed_area <- signed_area + (x_i * y_next - x_next * y_i)
  }
  
  # 调整顺时针方向（当面积为正时反转排序）
  if (signed_area > 0) {
    ordered_df <- ordered_df[c(1, n:2), ]
  }
  
  # 清理临时列并重置索引
  ordered_df$polar_angle <- NULL
  rownames(ordered_df) <- NULL
  
  return(ordered_df)
}


st_voronoi_point <- function(points){
  ## points must be POINT geometry
  # check for point geometry and execute if true
  if(!all(st_geometry_type(points) == "POINT")){
    stop("Input not  POINT geometries")
  }
  g = st_combine(st_geometry(points)) # make multipoint
  v = st_voronoi(g)
  v = st_collection_extract(v)
  return(v[unlist(st_intersects(points, v))])
}



correct_to_rectangle <- function(df) {
    
  # 定义原始坐标
  A <- unlist(df[1, c("x", "y")])
  B <- unlist(df[2, c("x", "y")])
  C <- unlist(df[3, c("x", "y")])
  D <- unlist(df[4, c("x", "y")])
  
  # 定义优化目标函数
  objective <- function(shift) {
    sum(shift^2)  # 最小化总平移量
  }
  
  # 定义非线性约束
  constraint <- function(shift) {
    # 解包位移向量（dxA, dyA, dxB, dyB, dxC, dyC, dxD, dyD）
    shifts <- matrix(shift, ncol = 2, byrow = TRUE)
    
    # 计算新坐标
    A_new <- A + shifts[1, ]
    B_new <- B + shifts[2, ]
    C_new <- C + shifts[3, ]
    D_new <- D + shifts[4, ]
    
    # 矩形约束条件
    vec_AB <- B_new - A_new
    vec_AD <- D_new - A_new
    vec_BC <- C_new - B_new
    
    # 正交性约束
    ortho <- sum(vec_AB * vec_AD)
    
    # 长度相等约束
    len_AB <- sqrt(sum(vec_AB^2))
    len_AD <- sqrt(sum(vec_AD^2))
    len_BC <- sqrt(sum(vec_BC^2))
    
    return(c(
      ortho,          # 正交性约束 (应为0)
      len_AB - len_AD, # 邻边长度相等 (应为0)
      len_BC - len_AD  # 对边长度相等 (应为0)
    ))
  }
  
  # 优化参数设置
  opts <- list(
    "algorithm" = "NLOPT_LN_AUGLAG",
    "xtol_rel" = 1e-6,
    "maxeval" = 1000,
    "local_opts" = list("algorithm" = "NLOPT_LN_COBYLA", "xtol_rel" = 1e-6)
  )
  
  # 初始猜测（零位移）
  init_shift <- rep(0, 8)
  
  # 运行优化
  solution <- nloptr(
    x0 = init_shift,
    eval_f = objective,
    eval_g_ineq = NULL,
    eval_g_eq = constraint,
    opts = opts,
    lb = rep(-Inf, 8),
    ub = rep(Inf, 8)
  )
  
  # 提取最优解
  optimal_shifts <- matrix(solution$solution, ncol = 2, byrow = TRUE)
  
  # 生成结果数据框
  result <- data.frame(
    orig_x = df$x,
    orig_y = df$y,
    new_x = df$x + optimal_shifts[,1],
    new_y = df$y + optimal_shifts[,2]
  )
  
  return(result)
}

transform_to_rectangle <- function(points) {
  # 检查输入是否为4个点
  if (nrow(points) != 4 || ncol(points) != 2) {
    stop("输入必须为4x2矩阵")
  }
  
  # 第一步：计算重心
  centroid <- colMeans(points)
  
  # 第二步：计算相对于重心的角度并顺时针排序
  angles <- atan2(points[, 2] - centroid[2], points[, 1] - centroid[1])
  order_idx <- order(angles, decreasing = TRUE)
  points_ordered <- points[order_idx, ]
  
  # 第三步：计算四条边的斜率
  slopes <- numeric(4)
  for (i in 1:4) {
    p1 <- points_ordered[i, ]
    p2 <- points_ordered[if (i == 4) 1 else i + 1, ]
    dx <- p2[1] - p1[1]
    dy <- p2[2] - p1[2]
    slopes[i] <- if (dx == 0) Inf else dy / dx
  }
  
  # 找到斜率绝对值最小的边
  min_slope_idx <- which.min(abs(unlist(slopes)))
  
  # 重新排序，使该边成为p1-p2
  if (min_slope_idx == 1) {
    points_reordered <- points_ordered
  } else if (min_slope_idx == 2) {
    points_reordered <- points_ordered[c(2, 3, 4, 1), ]
  } else if (min_slope_idx == 3) {
    points_reordered <- points_ordered[c(3, 4, 1, 2), ]
  } else {
    points_reordered <- points_ordered[c(4, 1, 2, 3), ]
  }
  
  # 定义旋转函数
  rotate <- function(point, center, angle) {
    dx <- point[1] - center[1]
    dy <- point[2] - center[2]
    new_x <- center[1] + dx * cos(angle) - dy * sin(angle)
    new_y <- center[2] + dx * sin(angle) + dy * cos(angle)
    c(new_x, new_y)
  }
  
  # 计算旋转角度，使p1-p2水平
  p1 <- unlist(points_reordered[1, ])
  p2 <- unlist(points_reordered[2, ])
  theta <- atan2(p2[2] - p1[2], p2[1] - p1[1])
  phi <- -theta
  
  # 旋转所有点
  points_rotated <- t(sapply(1:4, function(i) rotate(points_reordered[i, ], p1, phi)))
  
  # 获取旋转后的点
  p1_rot <- unlist(points_rotated[1, ])
  p2_rot <- unlist(points_rotated[2, ])
  p3_rot <- unlist(points_rotated[3, ])
  p4_rot <- unlist(points_rotated[4, ])
  
  # 第四步：调整p3和p4形成矩形
  y1 <- (p3_rot[2] + p4_rot[2]) / 2  # p3和p4的y坐标平均值
  p3_final <- c(p2_rot[1], y1)       # p3的x与p2对齐
  p4_final <- c(p1_rot[1], y1)       # p4的x与p1对齐
  
  # 返回最终的矩形点
  points_final <- rbind(p1_rot, p2_rot, p3_final, p4_final)
  return(points_final)
}


getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}

mirror_point <- function(A, B, P) {
  # 验证输入格式
  if (!all(sapply(list(A,B,P), length) == 2) || 
      !all(sapply(list(A,B,P), is.numeric))) {
    stop("输入必须为数值向量，格式：c(x,y)")
  }
  
  # 计算向量AB和AP
  AB <- B - A
  AP <- P - A
  
  # 检查直线有效性
  if (all(AB == 0)) {
    stop("点A和B重合，无法构成直线")
  }
  
  # 计算投影参数
  t <- sum(AB * AP) / sum(AB^2)
  
  # 计算投影点
  projection <- A + t * AB
  
  # 计算对称点
  Q <- 2 * projection - P
  
  return(round(Q, 2))  # 保留6位小数消除浮点误差
}

# 辅助函数：计算两向量的夹角（弧度）
angle_between_vectors <- function(v1, v2) {
  dot_product <- sum(v1 * v2)
  norm_v1 <- sqrt(sum(v1^2))
  norm_v2 <- sqrt(sum(v2^2))
  cos_theta <- dot_product / (norm_v1 * norm_v2)
  cos_theta <- max(min(cos_theta, 1), -1)  # 防止浮点误差
  theta <- acos(cos_theta)
  return(theta)
}

# 辅助函数：判断内角是否为凹角
is_concave_angle <- function(p1, p2, p3) {
  v1 <- p1 - p2
  v2 <- p3 - p2
  theta <- angle_between_vectors(v1, v2)
  # 计算叉积以判断转向
  cross_product <- v1[1] * v2[2] - v1[2] * v2[1]
  if (cross_product < 0) {
    theta <- 2 * pi - theta
  }
  return(theta > pi)
}

# 辅助函数：计算两条直线的交点
line_intersection <- function(line1, line2) {
  p1 <- line1[1, ]
  p2 <- line1[2, ]
  p3 <- line2[1, ]
  p4 <- line2[2, ]
  
  denom <- (p4[2] - p3[2]) * (p2[1] - p1[1]) - (p4[1] - p3[1]) * (p2[2] - p1[2])
  if (denom == 0) {
    stop("Lines are parallel and do not intersect.")
  }
  
  ua <- ((p4[1] - p3[1]) * (p1[2] - p3[2]) - (p4[2] - p3[2]) * (p1[1] - p3[1])) / denom
  intersection <- unlist(p1) + unlist(ua) * unlist(p2 - p1)
  return(intersection)
}

# 主函数：将凹四边形转换为凸四边形
convert_to_convex_quadrilateral <- function(points) {
  if (nrow(points) != 4 || ncol(points) != 2) {
    stop("Input must be a 4x2 matrix of points.")
  }
  
  # 识别凹角
  concave_index <- NULL
  for (i in 1:4) {
    p1 <- points[i, ]
    p2 <- points[(i %% 4) + 1, ]
    p3 <- points[(i + 1) %% 4 + 1, ]
    if (is_concave_angle(p1, p2, p3)) {
      concave_index <- (i %% 4) + 1
      break
    }
  }
  
  # 凹角相邻的点
  prev_index <- if (concave_index == 1) 4 else concave_index - 1
  next_index <- if (concave_index == 4) 1 else concave_index + 1
  last_index <- setdiff(c(1,2,3,4),c(prev_index,next_index,concave_index))
  prev_point <- points[prev_index, ]
  concave_point <- points[concave_index, ]
  next_point <- points[next_index, ]
  last_point <- points[last_index, ]
  # 定义两条直线并计算交点
  line1 <- rbind(prev_point, concave_point)
  line2 <- rbind(last_point, next_point)
  intersection_next <- line_intersection(line1, line2)
  line3 <- rbind(next_point, concave_point)
  line4 <- rbind(last_point, prev_point)
  intersection_prev <- line_intersection(line3, line4)
  # 形成新四边形
  new_points <- points
  new_points[prev_index, ] <- intersection_prev
  new_points[next_index, ] <- intersection_next
  
  return(new_points)
}


compare_strings <- function(str1, str2) {
  
  # 将字符串拆分为单个字符
  chars1 <- strsplit(str1, "")[[1]]
  chars2 <- strsplit(str2, "")[[1]]
  
  # 检查差异位置
  if (chars1[1] != chars2[1]) {
    return(TRUE)
  } else if (chars1[2] != chars2[2]) {
    return(FALSE)
  }
}

assign_direction <- function(df,row_index, medianx,mediany, pastlabel, condition) {
  # 验证输入条件
  if (!condition %in% c(TRUE, FALSE)) {
    stop("条件参数必须是'T'或'F'")
  }
  # 根据条件进行比较
  if (condition) {
    # 比较x值
    if (df[row_index, 1] > medianx) {
      df[row_index, 7] <- paste("R", strsplit(pastlabel, "")[[1]][2], sep = "")
    } else {
      df[row_index, 7] <- paste("L", strsplit(pastlabel, "")[[1]][2], sep = "")
    }
  } else {
    # 比较y值
    if (df[row_index, 2] > mediany) {
      df[row_index, 7] <- paste(strsplit(pastlabel, "")[[1]][1], "U", sep = "")
    } else {
      df[row_index, 7] <- paste(strsplit(pastlabel, "")[[1]][1], "D", sep = "")
    }
  }
  
  return(df)
}


which.median <- function(x){
  which.min(abs(x - median(x)))
} 

find_concave_vertex <- function(points) {
  # points: 一个4行2列的数据框或矩阵，每行是一个顶点的坐标，按顺时针顺序排列
  # 列顺序为x坐标和y坐标
  
  # 确保输入为4个点
  if (nrow(points) != 4) {
    stop("输入必须包含4个顶点")
  }
  
  # 为每个点计算叉积
  cross_products <- numeric(4)
  
  for (i in 1:4) {
    # 获取当前点、前一个点和后一个点
    prev_index <- ifelse(i == 1, 4, i - 1)
    next_index <- ifelse(i == 4, 1, i + 1)
    
    A <- points[prev_index, ]
    B <- points[i, ]
    C <- points[next_index, ]
    
    # 计算向量AB = (dx1, dy1) 和 BC = (dx2, dy2)
    dx1 <- B[1] - A[1]
    dy1 <- B[2] - A[2]
    dx2 <- C[1] - B[1]
    dy2 <- C[2] - B[2]
    
    # 计算叉积: AB x BC = (dx1 * dy2 - dy1 * dx2)
    cross_products[i] <- dx1 * dy2 - dy1 * dx2
  }
  
  # 在凹四边形中，凹点处的叉积应为正数，其余点为负数
  concave_index <- which(cross_products > 0)
  
  # 检查是否找到凹点
  if (length(concave_index) == 0) {
    stop("未找到凹点，请检查输入是否为凹四边形且顶点按顺时针排序")
  } else if (length(concave_index) > 1) {
    stop("找到多个凹点，输入可能不是简单凹四边形")
  }
  
  return(concave_index)
}

label_confirm <- function(df,rep1label1,rep1label2,rep2label1,rep2label2,meanx, meany) {
  pastelab1=paste(rep1label1,rep1label2, sep = "")
  pastelab2=paste(rep2label1,rep2label2, sep = "")
  if(grepl("LU",pastelab2)){
    df=rbind(df[2,],df[1,])
  }
  
  # 根据条件进行比较
  if (pastelab1=="LURU") {
    # 比较x值
    if (df[1, 1] > meanx) {
      df[1,7]="RU"
      if(df[2, 1] > meanx){
        df[2,7]="RD" 
      }else{
        df[2,7]="LD"
      } 
    }else{
      df[1,7]="LU"
      if(df[2, 1] > meanx){
        df[2,7]="RD" 
      }else{
        df[2,7]="LD"
      } 
    }
  }
  if (pastelab1=="LULD") {
    # 比较x值
    if (df[1, 2] > meany) {
      df[1,7]="LU"
      if(df[2, 2] > meany){
        df[2,7]="RU" 
      }else{
        df[2,7]="RD"
      } 
    }else{
      df[1,7]="LD"
      if(df[2, 2] > meany){
        df[2,7]="RU" 
      }else{
        df[2,7]="RD"
      } 
    }
  }  
  
  return(df)
}

find_max_triangle <- function(points) {
  # 检查输入是否为4x2矩阵
  if (nrow(points) != 4 || ncol(points) != 2) {
    stop("输入必须为4x2矩阵")
  }
  
  # 定义计算三角形面积的函数（鞋带公式）
  triangle_area <- function(p1, p2, p3) {
    area <- abs(0.5 * (
      p1[1] * (p2[2] - p3[2]) +
        p2[1] * (p3[2] - p1[2]) +
        p3[1] * (p1[2] - p2[2])
    ))
    return(area)
  }
  
  # 所有可能的三点组合和对应的剩余点
  combos <- list(
    list(triangle = c(1, 2, 3), remaining = 4),
    list(triangle = c(1, 2, 4), remaining = 3),
    list(triangle = c(1, 3, 4), remaining = 2),
    list(triangle = c(2, 3, 4), remaining = 1)
  )
  
  # 计算每个组合的面积
  max_area <- -Inf
  max_triangle_indices <- NULL
  remaining_index <- NULL
  
  for (combo in combos) {
    idx <- combo$triangle
    area <- triangle_area(points[idx[1], ], points[idx[2], ], points[idx[3], ])
    if (area > max_area) {
      max_area <- area
      max_triangle_indices <- idx
      remaining_index <- combo$remaining
    }
  }
  
  # 返回结果
  return(remaining_index = remaining_index)
}


detect_concave_quad <- function(points) {
  # 检查输入是否为4x2矩阵
  if (nrow(points) != 4 || ncol(points) != 2) {
    stop("输入必须为4x2矩阵")
  }
  
  points_ordered <- points
  indices_ordered <- c(1,2,3,4)
  
  # 线段相交检测函数
  line_intersect <- function(p1, p2, p3, p4) {
    # 计算线段p1-p2和p3-p4是否相交
    # 使用方向和跨立测试
    o1 <- function(x1, y1, x2, y2, x3, y3) {
      (y1 - y2) * (x3 - x2) - (x1 - x2) * (y3 - y2)
    }
    
    d1 <- o1(p1[1], p1[2], p2[1], p2[2], p3[1], p3[2])
    d2 <- o1(p1[1], p1[2], p2[1], p2[2], p4[1], p4[2])
    d3 <- o1(p3[1], p3[2], p4[1], p4[2], p1[1], p1[2])
    d4 <- o1(p3[1], p3[2], p4[1], p4[2], p2[1], p2[2])
    
    # 判断是否相交（排除端点相连的情况）
    return(d1 * d2 < 0 && d3 * d4 < 0)
  }
  
  # 检查四条边的相交情况
  is_concave <- FALSE
  concave_point <- NA
  
  # 定义四条边：1-2, 2-3, 3-4, 4-1
  edges <- list(
    c(1, 2, 3, 4), # 边1-2 vs 边3-4
    c(2, 3, 4, 1), # 边2-3 vs 边4-1
    c(3, 4, 1, 2), # 边3-4 vs 边1-2
    c(4, 1, 2, 3)  # 边4-1 vs 边2-3
  )
  
  for (i in 1:4) {
    p1 <- points_ordered[edges[[i]][1], ]
    p2 <- points_ordered[edges[[i]][2], ]
    p3 <- points_ordered[edges[[i]][3], ]
    p4 <- points_ordered[edges[[i]][4], ]
    
    if (line_intersect(p1, p2, p3, p4)) {
      is_concave <- TRUE
      break
    }
  }
  
  return(is_concave = is_concave)
}



######################################################################################################
######################################################################################################
######################################################################################################


library(terra)
library(sp)
library(rgeos)
library(sf)
library(tidyverse)
library(tigris)
library(lwgeom)
library(concaveman)
library(units)
library(LearnGeom)
library(mapview)
library(dplyr)
library(purrr)
library(tibble)
library(pracma)
library(ReinforcementLearning)
library(combinat)
library(nloptr)
library(sqldf)

a1=read.csv(file.choose())
data1=as.data.frame(a1)
colnames(data1)=c("x","y","no")
buff=2
data2=data1[data1$x>buff&data1$x<28,]
data2=data2[data2$y>buff&data2$y<18,]


#data3 is used for restore the result of tree-based UAI
data3=data1
buffer=rbind(c(2, 2), c(28, 2), c(28, 18), c(2, 18), c(2, 2)) %>% 
  list %>% 
  st_polygon %>% 
  st_sfc
roi=rbind(c(0.0, 0.0), c(30.0, 0.0), c(30.0, 20.0), c(0.0, 20.0), c(0.0, 0.0)) %>% 
  list %>% 
  st_polygon %>% 
  st_sfc

roisf=st_sf(roi)
combinall = data.frame()
cb = data.frame()
Polyunionfin=list()
tempx=c()
tempy=c()

####Define 3 sets,upper left(LU), lower left (LD), upper right (RU), and lower right (RD)
group1[1,1]="LD"
group1[1,2]="RD"
group1[2,1]="LU"
group1[2,2]="RD"
group1[3,1]="LU"
group1[3,2]="RU"
group1[4,1]="LD"
group1[4,2]="RU"

group2[1,1]="RU"
group2[1,2]="RD"
group2[2,1]="LU"
group2[2,2]="RD"
group2[3,1]="LD"
group2[3,2]="LU"
group2[4,1]="LD"
group2[4,2]="RU"

group3[1,1]="LU"
group3[1,2]="RU"
group3[2,1]="RU"
group3[2,2]="RD"
group3[3,1]="LD"
group3[3,2]="RD"
group3[4,1]="LD"
group3[4,2]="LU"

##Caclulate tree-based UAI
for (i in 1:length(data3$x)){
  x0=data3[i,]$x;y0=data3[i,]$y
  Near5=Near.f(x0,y0,data3)[1:5,]
  Near4=Near5[2:5,]
  data31$W[i]=w.f(x0,y0,Near4$x,Near4$y) #function 4
}
data3=data3[data3$x>buff&data3$x<28,]
data3=data3[data3$y>buff&data3$y<18,]
##Resule of tree-based UAI
mean(data3$W)


######Dot plot.(Due to the large amount of calculation, the dot plot takes 12 hours to complete)
x=seq(from=2,to=28,by=0.01)
y=seq(from=2,to=18,by=0.01)

for (i in 1:4164201) {
  k=i%%2601
  if(k==0){
    k=2601
  }
  tempx[i]=x[k]
}
  

for (i in 1:4164201) {
  k=floor((i-1)%/%2601)+1
  tempy[i]=y[k]
}
data6=tibble(x=tempx,y=tempy)
data6$no=1

data6=dplyr::setdiff(data6[,-3],data1[,-3])
data6$no=1

#timestart<-Sys.time()
#for(i in 1:4164092){
#  datal111 = rbind(data6[i,],data1)#combine trees in buffer area for UAI calculation
#  x0=data6[i,]$x;y0=data6[i,]$y
#  Near5=Near.f(x0,y0,datal111)[1:5,]
#  Near4=Near5[2:5,]
#  Near4=Near4[order(Near4$no),]
#  neb=paste(Near4$no[1],Near4$no[2],Near4$no[3],Near4$no[4],sep='')
#  data6[i,5]=w.f(x0,y0,Near4$x,Near4$y)#calculate UAI 
#  if (neb %in% k){
#    next
#  }
#  k[i]=neb
  #
#}   
#k <- k[!is.na(k)]
#timeend<-Sys.time()
#runningtime<-timeend-timestart
#print(runningtime)


###################Automated computing point-based UAI##########################

###Automatically generation of the RNN,this process takes about 30 minutes
timestart = Sys.time()
nc = st_read(system.file("shape/nc.shp", package="sf"), quiet = TRUE)[,15]
x=seq(from=2,to=28,by=0.1)
y=seq(from=2,to=18,by=0.1)
tempx=c()
tempy=c()
for (i in 1:(length(x)*length(y))) {
  k=i%%length(x)
  if(k==0){
    k=length(x)
  }
  tempx[i]=x[k]
}
datadot=data.frame(tempx)
colnames(datadot)=c("x")
for (i in 1:(length(x)*length(y))) {
  k=floor((i-1)%/%length(x))+1
  tempy[i]=y[k]
}
datadot$y=tempy
datadot1=dplyr::setdiff(datadot,data1[,-3])
datadot=dplyr::intersect(datadot,datadot1)
datadot$no=1

for(u in 1:100){
  neb=c()
  k=c()
  zoi=nc[1,]
  zoi=st_set_crs(zoi, NA)
  for(i in 1:length(datadot$x)){
    datal111 = rbind(datadot[i,],data1)#combine trees in buffer area for UAI calculation
    x0=datadot[i,]$x;y0=datadot[i,]$y
    Near5=Near.f(x0,y0,datal111)[1:5,]
    Near4=Near5[2:5,]
    Near4=Near4[order(Near4$no),]
    neb=paste(Near4$no[1],Near4$no[2],Near4$no[3],Near4$no[4],sep='')
    if (neb %in% k){
      next
    }
    k[i]=neb
  }   
  k <- k[!is.na(k)]
  f=substr(k,1,4)
  s=substr(k,5,8)
  t=substr(k,9,12)
  forth=substr(k,13,16)
  
  
  flist = as.data.frame(f)
  flist$s=s
  flist$t=t
  flist$forth=forth
  itimes=length(flist$f)
  
  
  for (o in 1:itimes) {
    for (z in 1:4) {
      select1=data1[which(data1$no %in% flist[o,z]),]
      elsepot = setdiff(unlist(flist[o,]),unlist(flist[o,z]))
      select11=data1[-which(data1$no %in% elsepot),]
      pot1=st_point(c(select1$x,select1$y),dim = "XY")
      pot1=SpatialPoints(coords = select1[,-3])
      pot11=SpatialPoints(coords = cbind(select11$x,select11$y))
      pot1 = st_as_sf(pot1)
      pot11 = st_as_sf(pot11)
      vor1 = st_voronoi_point(pot11)
      vor1 = st_as_sf(vor1)
      kkk=st_intersects(pot1, vor1, sparse = T) %>% 
        as.numeric()
      zone= vor1[kkk, ]
      if(z==1){
        zoiTemp=zone
      }
      zoiTemp=st_intersection(zone, zoiTemp)
    }
    zoi[o,]=zoiTemp
  }
  zoi=st_intersection(zoi,buffer)
  if(u==1){
    flisttemp=flist
    zoirecod=zoi
    finely_divided=zoi
  }
  if('POLYGON' %in% st_geometry_type(st_sym_difference(st_union(zoirecod,finely_divided),buffer), by_geometry = T)|
     'MULTIPOLYGON' %in% st_geometry_type(st_sym_difference(st_union(zoirecod,finely_divided),buffer), by_geometry = T)){
    otmit=st_sym_difference(st_union(zoirecod),buffer)
    otmit=st_cast(otmit, "MULTIPOLYGON") %>% st_cast("POLYGON")
    area <- st_area(otmit)
    if(length(which(area>0.00001))!=0){
      otmit=otmit[which(area>0.00001)]
      finely_divided=otmit[which(area<0.00001)]
      addpoints=st_sample(otmit, length(otmit), type="random")
      datadot=do.call(rbind, st_geometry(addpoints)) %>% 
        as_tibble()%>% setNames(c("x","y"))
      datadot$no=1
      if(u>1){
        flisttemp=rbind(flist,flisttemp)
        zoirecod=rbind(zoi,zoirecod)
      }
    }else{
      break
    }
  }else{
    break
  }
}
zoi=unique(zoirecod)
flist=unique(flisttemp)
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)







###Automatically Caclulate tree-based UAI,this process takes about 1 hour
timestart = Sys.time()
for(i in 1:length(flist[,1])){
####Determination of the four vertices of a quadrilateral
  local[1,] = data1[which(data1$no %in% flist[i,1]),]
  local[2,] = data1[which(data1$no %in% flist[i,2]),]
  local[3,] = data1[which(data1$no %in% flist[i,3]),]
  local[4,] = data1[which(data1$no %in% flist[i,4]),]
  
  
  x_median <- median(local[,1])
  y_median <- median(local[,2])
  local = sort_points_clockwise(local[,1:2])
  convex_check=local[,1:2]
  convex_check_result=is_convex_quadrilateral(convex_check)
  checkspecical=FALSE
  change=TRUE
  lookupvec=c("LU","RU","LD","RD")
  if(convex_check_result==T){
    x_median <- median(local[,1])
    y_median <- median(local[,2])
    local$LUindex=(local[,2]-y_median)-(local[,1]-x_median)
    local$RUindex=(local[,2]-y_median)+(local[,1]-x_median)
    local$LDindex=-(local[,2]-y_median)-(local[,1]-x_median)
    local$RDindex=(local[,1]-x_median)-(local[,2]-y_median)
    local[which.max(local$LUindex),7]="LU"
    local[which.max(local$LDindex),7]="LD"
    local[which.max(local$RUindex),7]="RU"
    local[which.max(local$RDindex),7]="RD"
    maxindex=c(which.max(local$LUindex),which.max(local$RUindex),which.max(local$LDindex),which.max(local$RDindex))
    if(length(unique(maxindex))!=4){
      if(length(unique(maxindex))==3){
        leck=setdiff(c(1,2,3,4),maxindex)
        repect=getmode(maxindex)
        repindex=which(maxindex==repect)
        reptype=lookupvec[repindex]
        if(is.na(sum(match(c("LU","RU"),reptype)))==F){
          if(local[repect,1]>=local[leck,1]){
            local[repect,7]="RU"
            local[leck,7]="LU"
            local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
            convex_check_result=is_convex_quadrilateral(local[,1:2])
          }else{
            local[repect,7]="LU"
            local[leck,7]="RU"
            local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
            convex_check_result=is_convex_quadrilateral(local)
          }
        }
        if(is.na(sum(match(c("LD","RD"),reptype)))==F){
          if(local[repect,1]>=local[leck,1]){
            local[repect,7]="RD"
            local[leck,7]="LD"
            local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
            convex_check_result=is_convex_quadrilateral(local)
          }else{
            local[repect,7]="LD"
            local[leck,7]="RD"
            local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
            convex_check_result=is_convex_quadrilateral(local)
          }
        }
        if(is.na(sum(match(c("LU","LD"),reptype)))==F){
          if(local[repect,2]>=local[leck,2]){
            local[repect,7]="LU"
            local[leck,7]="LD"
            local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
            convex_check_result=is_convex_quadrilateral(local)
          }else{
            local[repect,7]="LD"
            local[leck,7]="LU"
            local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
            convex_check_result=is_convex_quadrilateral(local)
          }
        }
        if(is.na(sum(match(c("RU","RD"),reptype)))==F){
          if(local[repect,2]>=local[leck,2]){
            local[repect,7]="RU"
            local[leck,7]="RD"
            local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
            convex_check_result=is_convex_quadrilateral(local)
          }else{
            local[repect,7]="RD"
            local[leck,7]="RU"
            local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
            convex_check_result=is_convex_quadrilateral(local)
          }
        }
      }
      if(length(unique(maxindex))==2){
        repect=unique(maxindex)
        rep1label1=lookupvec[which(maxindex == repect[1])[1]] 
        rep1label2=lookupvec[which(maxindex == repect[1])[2]] 
        rep2label1=lookupvec[which(maxindex == repect[2])[1]] 
        rep2label2=lookupvec[which(maxindex == repect[2])[2]]
        leck=setdiff(c(1,2,3,4),maxindex)
        hull_points=rbind(local[repect[1],],local[repect[2],])
        hull_points_x_median=median(local[,1])
        hull_points_y_median=median(local[,2])
        hull_points=label_confirm(hull_points,rep1label1,rep1label2,rep2label1,rep2label2,hull_points_x_median, hull_points_y_median)
        currentlaber=hull_points[,7]
        stillleck=setdiff(c(1,2,3,4),match(currentlaber,lookupvec))
        lecklaber=lookupvec[stillleck]
        concave_point=rbind(local[leck[1],],local[leck[2],])
        concave_point_x_mean=mean(concave_point[,1])
        concave_point_y_mean=mean(concave_point[,2])
        
        if((grepl("R",lecklaber[1])|grepl("R",lecklaber[2]))&(grepl("L",lecklaber[1])|grepl("L",lecklaber[2]))){
          if((grepl("U",lecklaber[1])|grepl("U",lecklaber[2]))&(grepl("D",lecklaber[1])|grepl("D",lecklaber[2]))){
            lecklaberpas=paste(lecklaber[1],lecklaber[2],sep = "")
            if(lecklaberpas=="RULD"|lecklaberpas=="LDRU"){
              if(max(rank(concave_point$RUindex))!=2){
                if(concave_point[1,1]>concave_point[2,1]){
                  concave_point[1,7]="RU"
                  concave_point[2,7]="LD"
                }else{
                  concave_point[1,7]="LD"
                  concave_point[2,7]="RU"
                }
              }else{
                concave_point[which.max(concave_point$RUindex),7]="RU"
                concave_point[which.min(concave_point$RUindex),7]="LD" 
              }
            }
            if(lecklaberpas=="RDLU"|lecklaberpas=="LURD"){
              if(max(rank(concave_point$LUindex))!=2){
                if(concave_point[1,1]>concave_point[2,1]){
                  concave_point[1,7]="RD"
                  concave_point[2,7]="LU"
                }else{
                  concave_point[1,7]="LU"
                  concave_point[2,7]="RD"
                }
              }else{
                concave_point[which.max(concave_point$LUindex),7]="LU"
                concave_point[which.min(concave_point$LUindex),7]="RD"
              }
            }
            
          }else{
            if(concave_point[1,1]!=concave_point[2,1]){
              if(concave_point_x_mean>concave_point[1,1]){
                concave_point[2,7]=lecklaber[grep("R",lecklaber)]
                concave_point[1,7]=lecklaber[grep("L",lecklaber)]
              }
              if(concave_point_x_mean<concave_point[1,1]){
                concave_point[1,7]=lecklaber[grep("R",lecklaber)]
                concave_point[2,7]=lecklaber[grep("L",lecklaber)]
              }
            }else{
              if(sum(grepl("U",lecklaber[1]), na.rm = TRUE)==1&sum(grepl("U",lecklaber[2]), na.rm = TRUE)==1){ 
                if(hull_points_y_mean>concave_point[1,2]){
                  hull_points[3,7]=lecklaber[grep("R",lecklaber)]
                  concave_point[1,7]=lecklaber[grep("L",lecklaber)]
                }
                if(hull_points_y_mean<concave_point[1,2]){
                  hull_points[3,7]=lecklaber[grep("L",lecklaber)]
                  concave_point[1,7]=lecklaber[grep("R",lecklaber)]
                }  
              }
              else{
                if(hull_points_y_mean>concave_point[1,2]){
                  concave_point[1,7]=lecklaber[grep("R",lecklaber)]
                  concave_point[2,7]=lecklaber[grep("L",lecklaber)]
                }
                if(hull_points_y_mean<concave_point[1,2]){
                  concave_point[1,7]=lecklaber[grep("L",lecklaber)]
                  concave_point[2,7]=lecklaber[grep("R",lecklaber)]
                }  
              }
            } 
          }
        }else{
          if((grepl("R",lecklaber[1])&grepl("R",lecklaber[2]))|(grepl("L",lecklaber[1])&grepl("L",lecklaber[2]))){
            if(concave_point_y_mean>concave_point[2,2]){
              concave_point[1,7]=lecklaber[grep("U",lecklaber)]
              concave_point[2,7]=lecklaber[grep("D",lecklaber)]
            }
            if(concave_point_y_mean<concave_point[2,2]){
              concave_point[1,7]=lecklaber[grep("D",lecklaber)]
              concave_point[2,7]=lecklaber[grep("U",lecklaber)]
            }
            if(concave_point[2,2]==concave_point[1,2]){
              if(concave_point[2,1]<concave_point[1,1]){
                concave_point[2,7]=lecklaber[grep("U",lecklaber)]
                concave_point[1,7]=lecklaber[grep("D",lecklaber)]
              }
              else{
                concave_point[1,7]=lecklaber[grep("U",lecklaber)]
                concave_point[2,7]=lecklaber[grep("D",lecklaber)]
              }
            }
          }
        }
        
        local=rbind(hull_points,concave_point)
        local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
        convex_check_result=is_convex_quadrilateral(local)
      }
    }
    if(local[grep("RU",local$V7),1]<=local[grep("LU",local$V7),1]|local[grep("RD",local$V7),1]<=local[grep("LD",local$V7),1]|
       local[grep("RU",local$V7),2]<=local[grep("RD",local$V7),2]|local[grep("LU",local$V7),2]<=local[grep("LD",local$V7),2]){
      if(detect_concave_quad(local[,1:2])==T){
        if(local[grep("RU",local$V7),1]<=local[grep("LU",local$V7),1]){
          local[grep("RU",local$V7),7]="ru"
          local[grep("LU",local$V7),7]="RU"
          local[grep("ru",local$V7),7]="LU"
          local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
          if(local[grep("RU",local$V7),1]<local[grep("LU",local$V7),1]){
            convex_check_result=is_convex_quadrilateral(local)
          }else{
            convex_check_result=F
            change=F
            ccave=find_max_triangle(local[,1:2])
            if(4-ccave==1){
              local=rbind(local[4,],local[1:3,])
            }
            if(4-ccave==2){
              local=rbind(local[3:4,],local[1:2,])
            }
            if(4-ccave==3){
              local=rbind(local[2:4,],local[1,])
            }
          }
        }
        if(local[grep("RD",local$V7),1]<=local[grep("LD",local$V7),1]){
          local[grep("RD",local$V7),7]="rd"
          local[grep("LD",local$V7),7]="RD"
          local[grep("rd",local$V7),7]="LD"
          local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
          if(local[grep("RD",local$V7),1]<local[grep("LD",local$V7),1]){
            convex_check_result=is_convex_quadrilateral(local)
          }else{
            convex_check_result=F
            change=F
            ccave=find_max_triangle(local[,1:2])
            if(4-ccave==1){
              local=rbind(local[4,],local[1:3,])
            }
            if(4-ccave==2){
              local=rbind(local[3:4,],local[1:2,])
            }
            if(4-ccave==3){
              local=rbind(local[2:4,],local[1,])
            }
          }
        }
        if(local[grep("RU",local$V7),2]<=local[grep("RD",local$V7),2]){
          local[grep("RU",local$V7),7]="ru"
          local[grep("RD",local$V7),7]="RU"
          local[grep("ru",local$V7),7]="RD"
          local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
          if(local[grep("RU",local$V7),2]<local[grep("RD",local$V7),2]){
            convex_check_result=is_convex_quadrilateral(local)
          }else{
            convex_check_result=F
            change=F
            ccave=find_max_triangle(local[,1:2])
            if(4-ccave==1){
              local=rbind(local[4,],local[1:3,])
            }
            if(4-ccave==2){
              local=rbind(local[3:4,],local[1:2,])
            }
            if(4-ccave==3){
              local=rbind(local[2:4,],local[1,])
            }
          }
        }
        if(local[grep("LU",local$V7),2]<=local[grep("LD",local$V7),2]){
          local[grep("LU",local$V7),7]="lu"
          local[grep("LD",local$V7),7]="LU"
          local[grep("lu",local$V7),7]="LD"
          local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
          if(local[grep("LU",local$V7),2]<local[grep("LD",local$V7),2]){
            convex_check_result=is_convex_quadrilateral(local)
          }else{
            convex_check_result=F
            change=F
            ccave=find_max_triangle(local[,1:2])
            if(4-ccave==1){
              local=rbind(local[4,],local[1:3,])
            }
            if(4-ccave==2){
              local=rbind(local[3:4,],local[1:2,])
            }
            if(4-ccave==3){
              local=rbind(local[2:4,],local[1,])
            }
          }
          
        }
      }else{
        checkspecical=TRUE
        convex_check_result=F
      }
    }else{
      if(convex_check_result==F){
        change=F
        dcquad=detect_concave_quad(local[,1:2])
        if(dcquad==T){
          if(min(rank(local[,1]))!=1){
            local[grep("LD",local$V7),7]="ru"
            local[grep("RU",local$V7),7]="LD"
            local[grep("ru",local$V7),7]="RU"
            local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
          }
        }else{
          ccave=find_max_triangle(local[,1:2])
          if(4-ccave==1){
            local=rbind(local[4,],local[1:3,])
          }
          if(4-ccave==2){
            local=rbind(local[3:4,],local[1:2,])
          }
          if(4-ccave==3){
            local=rbind(local[2:4,],local[1,])
          }
        }
      }
    }
  }
  if(change==T){
    if(convex_check_result==F){
      if(checkspecical==F){
        concave_indices <- find_max_triangle(local[,1:2])
        concave_point <- local[concave_indices, ]
        all_indices <- 1:4
        hull_index <- setdiff(all_indices, concave_indices)
        hull_points <- local[hull_index, ]
      }
      else{
        if(local[grep("RU",local$V7),1]<=local[grep("LU",local$V7),1]&local[grep("RU",local$V7),2]<=mean(local[,2])){
          concave_point <- local[grep("RU",local$V7), -7]
          hull_points <- local[-grep("RU",local$V7), -7]
        }
        if(local[grep("RU",local$V7),1]<=local[grep("LU",local$V7),1]&local[grep("RU",local$V7),2]>mean(local[,2])){
          concave_point <- local[grep("LU",local$V7), -7]
          hull_points <- local[-grep("LU",local$V7), -7]
        }
        if(local[grep("RD",local$V7),1]<=local[grep("LD",local$V7),1]&local[grep("RD",local$V7),2]<=mean(local[,2])){
          concave_point <- local[grep("RD",local$V7), -7]
          hull_points <- local[-grep("RD",local$V7), -7]
        }
        if(local[grep("RD",local$V7),1]<=local[grep("LD",local$V7),1]&local[grep("RD",local$V7),2]>mean(local[,2])){
          concave_point <- local[grep("LD",local$V7), -7]
          hull_points <- local[-grep("LD",local$V7), -7]
        }
        if(local[grep("RU",local$V7),2]<=local[grep("RD",local$V7),2]&local[grep("RU",local$V7),2]<=mean(local[,2])){
          concave_point <- local[grep("RU",local$V7), -7]
          hull_points <- local[-grep("RU",local$V7), -7]
        }
        if(local[grep("RU",local$V7),2]<=local[grep("RD",local$V7),2]&local[grep("RU",local$V7),2]>mean(local[,2])){
          concave_point <- local[grep("RD",local$V7), -7]
          hull_points <- local[-grep("RD",local$V7), -7]
        }
        if(local[grep("LU",local$V7),2]<=local[grep("LD",local$V7),2]&local[grep("LU",local$V7),2]<=mean(local[,2])){
          concave_point <- local[grep("LU",local$V7),-7 ]
          hull_points <- local[-grep("LU",local$V7),-7 ]
        }
        if(local[grep("LU",local$V7),2]<=local[grep("LD",local$V7),2]&local[grep("LU",local$V7),2]>mean(local[,2])){
          concave_point <- local[grep("LD",local$V7),-7 ]
          hull_points <- local[-grep("LD",local$V7),-7 ]
        }
      }
      hull_points_x_median <- median(hull_points[,1])
      hull_points_y_median <- median(hull_points[,2])
      hull_points$LUindex=(hull_points[,2]+hull_points_y_median)-(hull_points[,1]-hull_points_x_median)
      hull_points$RUindex=(hull_points[,2]-hull_points_y_median)+(hull_points[,1]-hull_points_x_median)
      hull_points$LDindex=-(hull_points[,2]-hull_points_y_median)-(hull_points[,1]-hull_points_x_median)
      hull_points$RDindex=(hull_points[,1]-hull_points_x_median)-(hull_points[,2]-hull_points_y_median)
      if(max(hull_points$RUindex)+max(hull_points$LDindex)<0.0000001|max(hull_points$LUindex)+max(hull_points$RDindex)<0.000001){
        if(max(hull_points$RUindex)+max(hull_points$LDindex)<0.0000001){
          if(concave_point[1,2]<=min(hull_points$y)){
            hull_points[which.min(hull_points$x),7]="LD"
            hull_points[which.median(hull_points$x),7]="LU"
            hull_points[which.max(hull_points$x),7]="RU"
            concave_point[1,7]="RD"
          }else{
            hull_points[which.min(hull_points$x),7]="LU"
            hull_points[which.median(hull_points$x),7]="LD"
            hull_points[which.max(hull_points$x),7]="RD"
            concave_point[1,7]="RU"
          }
        }else{
          if(concave_point[1,2]>=max(hull_points$y)){
            hull_points[which.min(hull_points$x),7]="LU"
            hull_points[which.median(hull_points$x),7]="LD"
            hull_points[which.max(hull_points$x),7]="RD"
            concave_point[1,7]="RU"
          }else{
            hull_points[which.min(hull_points$x),7]="LD"
            hull_points[which.median(hull_points$x),7]="LU"
            hull_points[which.max(hull_points$x),7]="RU"
            concave_point[1,7]="RD"
          }
        }
      }
      else{
        hull_points[which.max(hull_points$LUindex),7]="LU"
        hull_points[which.max(hull_points$LDindex),7]="LD"
        hull_points[which.max(hull_points$RUindex),7]="RU"
        hull_points[which.max(hull_points$RDindex),7]="RD"
        concave_point[,3:6] <- hull_points[1,3:6]#avoid error
        maxindex1=c(which.max(hull_points$LUindex),which.max(hull_points$RUindex),which.max(hull_points$LDindex),which.max(hull_points$RDindex))
        lookupvec=c("LU","RU","LD","RD")
        if(length(unique(maxindex1))==3){
          leck=setdiff(c(1,2,3,4),maxindex1)
          repect=getmode(maxindex1)
          repindex=which(maxindex1==repect)
          reptype=lookupvec[repindex]
          if(is.na(sum(match(c("LU","RU"),reptype)))==F){
            if(hull_points[repect,1]>=concave_point[1,1]){
              if(hull_points[repect,1]==concave_point[1,1]){
                if(hull_points[repect,1]>median(local[,1])){
                  hull_points[repect,7]="RU"
                  concave_point[1,7]="LU"
                  local=rbind(hull_points,concave_point)
                }else{
                  hull_points[repect,7]="LU"
                  concave_point[1,7]="RU"
                  local=rbind(hull_points,concave_point)
                }
              }else{
                hull_points[repect,7]="RU"
                concave_point[1,7]="LU"
                local=rbind(hull_points,concave_point)
                local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                convex_check_result=is_convex_quadrilateral(local)
                if(convex_check_result==F){
                  local=rbind(hull_points,concave_point)
                } 
              }
            }else{
              hull_points[repect,7]="LU"
              concave_point[1,7]="RU"
              local=rbind(hull_points,concave_point)
              local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
              convex_check_result=is_convex_quadrilateral(local)
              if(convex_check_result==F){
                local=rbind(hull_points,concave_point)
              }
            }
          }
          if(is.na(sum(match(c("LD","RD"),reptype)))==F){
            if(hull_points[repect,1]>=concave_point[1,1]){
              if(hull_points[repect,1]==concave_point[1,1]){
                if(hull_points[repect,1]>median(local[,1])){
                  hull_points[repect,7]="RD"
                  concave_point[1,7]="LD"
                  local=rbind(hull_points,concave_point)
                }else{
                  hull_points[repect,7]="LD"
                  concave_point[1,7]="RD"
                  local=rbind(hull_points,concave_point)
                }
              }else{
                hull_points[repect,7]="RD"
                concave_point[1,7]="LD"
                local=rbind(hull_points,concave_point)
                local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                convex_check_result=is_convex_quadrilateral(local)
                if(convex_check_result==F){
                  local=rbind(hull_points,concave_point)
                }
              }
            }else{
              hull_points[repect,7]="LD"
              concave_point[1,7]="RD"
              local=rbind(hull_points,concave_point)
              local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
              convex_check_result=is_convex_quadrilateral(local)
              if(convex_check_result==F){
                local=rbind(hull_points,concave_point)
              }
            }
          }
          if(is.na(sum(match(c("LU","LD"),reptype)))==F){
            if(hull_points[repect,2]>=concave_point[1,2]){
              if(hull_points[repect,2]==concave_point[1,2]){
                if(hull_points[repect,2]>median(local[,2])){
                  hull_points[repect,7]="LU"
                  concave_point[1,7]="LD"
                  local=rbind(hull_points,concave_point)
                }else{
                  hull_points[repect,7]="LD"
                  concave_point[1,7]="LU"
                  local=rbind(hull_points,concave_point)
                }
              }else{
                hull_points[repect,7]="LU"
                concave_point[1,7]="LD"
                local=rbind(hull_points,concave_point)
                local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                convex_check_result=is_convex_quadrilateral(local)
                if(convex_check_result==F){
                  local=rbind(hull_points,concave_point)
                }
              }
            }else{
              hull_points[repect,7]="LD"
              concave_point[1,7]="LU"
              local=rbind(hull_points,concave_point)
              local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
              convex_check_result=is_convex_quadrilateral(local)
              if(convex_check_result==F){
                local=rbind(hull_points,concave_point)
              }
            }
          }
          if(is.na(sum(match(c("RU","RD"),reptype)))==F){
            if(hull_points[repect,2]>=concave_point[1,2]){
              if(hull_points[repect,2]==concave_point[1,2]){
                if(hull_points[repect,2]>median(local[,2])){
                  hull_points[repect,7]="RU"
                  concave_point[1,7]="RD"
                  local=rbind(hull_points,concave_point)
                }else{
                  hull_points[repect,7]="RD"
                  concave_point[1,7]="RU"
                  local=rbind(hull_points,concave_point)
                }  
              }else{
                hull_points[repect,7]="RU"
                concave_point[1,7]="RD"
                local=rbind(hull_points,concave_point)
                local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                convex_check_result=is_convex_quadrilateral(local)
                if(convex_check_result==F){
                  local=rbind(hull_points,concave_point)
                } 
              }
            }else{
              hull_points[repect,7]="RD"
              concave_point[1,7]="RU"
              local=rbind(hull_points,concave_point)
              local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
              convex_check_result=is_convex_quadrilateral(local)
              if(convex_check_result==F){
                local=rbind(hull_points,concave_point)
              }
            }
          }
        }
        if(length(unique(maxindex1))==2){
          repect=unique(maxindex1)
          rep1label1=lookupvec[which(maxindex1 == repect[1])[1]] 
          rep1label2=lookupvec[which(maxindex1 == repect[1])[2]] 
          rep2label1=lookupvec[which(maxindex1 == repect[2])[1]] 
          rep2label2=lookupvec[which(maxindex1 == repect[2])[2]]
          rep1type=compare_strings(rep1label1,rep1label2)
          rep2type=compare_strings(rep2label1,rep2label2)
          leck=setdiff(c(1,2,3),maxindex1)
          hull_points=assign_direction(hull_points, repect[1], x_median,y_median, rep1label1, rep1type)
          hull_points=assign_direction(hull_points, repect[2], x_median,y_median, rep2label1, rep2type)
          currentlaber=hull_points[-leck,7]
          stillleck=setdiff(c(1,2,3,4),match(currentlaber,lookupvec))
          lecklaber=lookupvec[stillleck]
          if((grepl("R",lecklaber[1])|grepl("R",lecklaber[2]))&(grepl("L",lecklaber[1])|grepl("L",lecklaber[2]))){
            if(hull_points[leck,1]!=concave_point[1,1]){
              if(hull_points[leck,1]>concave_point[1,1]){
                hull_points[leck,7]=lecklaber[grep("R",lecklaber)]
                concave_point[1,7]=lecklaber[grep("L",lecklaber)]
                local=rbind(hull_points,concave_point)
                local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                convex_check_result=is_convex_quadrilateral(local)
                if(convex_check_result==F){
                  local=rbind(hull_points,concave_point)
                }
              }
              if(hull_points[leck,1]<concave_point[1,1]){
                hull_points[leck,7]=lecklaber[grep("L",lecklaber)]
                concave_point[1,7]=lecklaber[grep("R",lecklaber)]
                local=rbind(hull_points,concave_point)
                local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                convex_check_result=is_convex_quadrilateral(local)
                if(convex_check_result==F){
                  local=rbind(hull_points,concave_point)
                }
              }
            }else{
              if(sum(grepl("U",lecklaber), na.rm = TRUE)==1&sum(grepl("U",lecklaber), na.rm = TRUE)==1){ 
                if(hull_points[leck,2]>concave_point[1,2]){
                  hull_points[leck,7]=lecklaber[grep("U",lecklaber)]
                  concave_point[1,7]=lecklaber[grep("D",lecklaber)]
                  local=rbind(hull_points,concave_point)
                  local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                  convex_check_result=is_convex_quadrilateral(local)
                  if(convex_check_result==F){
                    local=rbind(hull_points,concave_point)
                  }
                }
                if(hull_points[leck,2]<concave_point[1,2]){
                  hull_points[leck,7]=lecklaber[grep("D",lecklaber)]
                  concave_point[1,7]=lecklaber[grep("U",lecklaber)]
                  local=rbind(hull_points,concave_point)
                  local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                  convex_check_result=is_convex_quadrilateral(local)
                  if(convex_check_result==F){
                    local=rbind(hull_points,concave_point)
                  }
                }  
              }else{
                if(hull_points[leck,2]>concave_point[1,2]){
                  hull_points[leck,7]=lecklaber[grep("L",lecklaber)]
                  concave_point[1,7]=lecklaber[grep("R",lecklaber)]
                  local=rbind(hull_points,concave_point)
                  local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                  convex_check_result=is_convex_quadrilateral(local)
                  if(convex_check_result==F){
                    local=rbind(hull_points,concave_point)
                  }
                }
                if(hull_points[leck,2]<concave_point[1,2]){
                  hull_points[leck,7]=lecklaber[grep("R",lecklaber)]
                  concave_point[1,7]=lecklaber[grep("L",lecklaber)]
                  local=rbind(hull_points,concave_point)
                  local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                  convex_check_result=is_convex_quadrilateral(local)
                  if(convex_check_result==F){
                    local=rbind(hull_points,concave_point)
                  }
                }  
              }
            }
          }
          else{
            if((grepl("R",lecklaber[1])&grepl("R",lecklaber[2]))|(grepl("L",lecklaber[1])&grepl("L",lecklaber[2]))){
              if(hull_points[leck,2]>concave_point[1,2]){
                hull_points[leck,7]=lecklaber[grep("U",lecklaber)]
                concave_point[1,7]=lecklaber[grep("D",lecklaber)]
                local=rbind(hull_points,concave_point)
                local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                convex_check_result=is_convex_quadrilateral(local)
                if(convex_check_result==F){
                  local=rbind(hull_points,concave_point)
                }
              }
              if(hull_points[leck,2]<concave_point[1,2]){
                hull_points[leck,7]=lecklaber[grep("D",lecklaber)]
                concave_point[1,7]=lecklaber[grep("U",lecklaber)]
                local=rbind(hull_points,concave_point)
                local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                convex_check_result=is_convex_quadrilateral(local)
                if(convex_check_result==F){
                  local=rbind(hull_points,concave_point)
                }
              }
              if(hull_points[leck,2]==concave_point[1,2]){
                if(hull_points[leck,1]<concave_point[1,1]){
                  hull_points[leck,7]=lecklaber[grep("U",lecklaber)]
                  concave_point[1,7]=lecklaber[grep("D",lecklaber)]
                  local=rbind(hull_points,concave_point)
                  local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                  convex_check_result=is_convex_quadrilateral(local)
                  if(convex_check_result==F){
                    local=rbind(hull_points,concave_point)
                  }
                }
                else{
                  hull_points[leck,7]=lecklaber[grep("D",lecklaber)]
                  concave_point[1,7]=lecklaber[grep("U",lecklaber)]
                  local=rbind(hull_points,concave_point)
                  local=rbind(local[grep("LU",local$V7),],local[grep("LD",local$V7),],local[grep("RD",local$V7),],local[grep("RU",local$V7),])
                  convex_check_result=is_convex_quadrilateral(local)
                  if(convex_check_result==F){
                    local=rbind(hull_points,concave_point)
                  }
                }
              }
            }
          }
        }
      }
    }
  }
  
  
  
  ###########Superposition of three sets of double arcs
  LU=local[which(local$V7=="LU"),]
  LD=local[which(local$V7=="LD"),]
  RU=local[which(local$V7=="RU"),]
  RD=local[which(local$V7=="RD"),]
  LULD=rbind(LU,LD)
  RURD=rbind(RU,RD)
  RULD=rbind(RU,LD)
  LURU=rbind(LU,RU)
  LURD=rbind(LU,RD)
  LDRD=rbind(LD,RD)
  
  ovrlopg1=list()
  ovrlopg2=list()
  ovrlopg4=list()
  
  for(z in 1:4){
    point1=local[which(local$V7==group1[z,1]),]
    point2=local[which(local$V7==group1[z,2]),]
    ovrlopg1[[z]] = circle(point1,point2)
  }
  overlop11= st_sfc(ovrlopg1)
  overlop11= st_sf(overlop11)
  grp1 = overlop11 %>% st_set_precision(1000) %>% st_intersection()
  ##plot(grp1["n.overlaps"])
  ##grp1 = st_intersection(overlop11)
  
  for(z in 1:4){
    point1=local[which(local$V7==group2[z,1]),]
    point2=local[which(local$V7==group2[z,2]),]
    ovrlopg2[[z]] = circle(point1,point2)
  }
  overlop22= st_sfc(ovrlopg2)
  overlop22= st_sf(overlop22)
  grp2 = overlop22 %>% st_set_precision(1000) %>% st_intersection()
  
  
  for(z in 1:4){
    point1=local[which(local$V7==group3[z,1]),]
    point2=local[which(local$V7==group3[z,2]),]
    ovrlopg4[[z]] = circle(point1,point2)
  }
  overlop44= st_sfc(ovrlopg4)
  overlop44= st_sf(overlop44)
  grp3 = overlop44 %>% st_set_precision(1000) %>% st_intersection()
  
  LULD_PS=extend_point(c(LU[1,1],LU[1,2]),c(LD[1,1],LD[1,2]))
  LDLU_PS=extend_point(c(LD[1,1],LD[1,2]),c(LU[1,1],LU[1,2]))
  LURU_PS=extend_point(c(LU[1,1],LU[1,2]),c(RU[1,1],RU[1,2]))
  RULU_PS=extend_point(c(RU[1,1],RU[1,2]),c(LU[1,1],LU[1,2]))
  LURD_PS=extend_point(c(LU[1,1],LU[1,2]),c(RD[1,1],RD[1,2]))
  RDLU_PS=extend_point(c(RD[1,1],RD[1,2]),c(LU[1,1],LU[1,2]))
  RURD_PS=extend_point(c(RU[1,1],RU[1,2]),c(RD[1,1],RD[1,2]))
  RDRU_PS=extend_point(c(RD[1,1],RD[1,2]),c(RU[1,1],RU[1,2]))
  RULD_PS=extend_point(c(RU[1,1],RU[1,2]),c(LD[1,1],LD[1,2]))
  LDRU_PS=extend_point(c(LD[1,1],LD[1,2]),c(RU[1,1],RU[1,2]))
  LDRD_PS=extend_point(c(LD[1,1],LD[1,2]),c(RD[1,1],RD[1,2]))
  RDLD_PS=extend_point(c(RD[1,1],RD[1,2]),c(LD[1,1],LD[1,2]))
  
  ################################Generation of spatial configuration areas and multifactor overlay
  if(convex_check_result==F){
    if(local[4,7]=="RU"){
      
      group1Poly1 <- CreatePolygon(c(LURD_PS[1],LURD_PS[2]), c(LU$x,LU$y), c(LURU_PS[1],LURU_PS[2]),c(LURD_PS[1],LURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group1Poly2 <- CreatePolygon(c(LDRD_PS[1],LDRD_PS[2]), c(LD$x,LD$y), c(LDRU_PS[1],LDRU_PS[2]),c(LDRD_PS[1],LDRD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group1Poly3 <- CreatePolygon(c(RULU_PS[1],RULU_PS[2]), c(RU$x,RU$y), c(RULD_PS[1],RULD_PS[2]),c(RDLD_PS[1],RDLD_PS[2]),c(RD$x,RD$y),c(RDLU_PS[1],RDLU_PS[2]),c(RULU_PS[1],RULU_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly1 <- CreatePolygon(c(LDRU_PS[1],LDRU_PS[2]), c(LD$x,LD$y), c(LDLU_PS[1],LDLU_PS[2]),c(LDRU_PS[1],LDRU_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly2 <- CreatePolygon(c(RDRU_PS[1],RDRU_PS[2]), c(RD$x,RD$y), c(RDLU_PS[1],RDLU_PS[2]),c(RDRU_PS[1],RDRU_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly3 <- CreatePolygon(c(RURD_PS[1],RURD_PS[2]), c(RU$x,RU$y), c(RULD_PS[1],RULD_PS[2]),c(LULD_PS[1],LULD_PS[2]),c(LU$x,LU$y),c(LURD_PS[1],LURD_PS[2]),c(RURD_PS[1],RURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly1 <- CreatePolygon(c(RDLD_PS[1],RDLD_PS[2]), c(RD$x,RD$y), c(RDRU_PS[1],RDRU_PS[2]),c(RDLD_PS[1],RDLD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly2 <- CreatePolygon(c(LURU_PS[1],LURU_PS[2]), c(LU$x,LU$y), c(LULD_PS[1],LULD_PS[2]),c(LURU_PS[1],LURU_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly3 <- CreatePolygon(c(RURD_PS[1],RURD_PS[2]), c(RU$x,RU$y), c(RULU_PS[1],RULU_PS[2]),c(LDLU_PS[1],LDLU_PS[2]),c(LD$x,LD$y),c(LDRD_PS[1],LDRD_PS[2]),c(RURD_PS[1],RURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      
      group1Poly1=sf::st_buffer(group1Poly1, dist = 0)
      group1Poly2=sf::st_buffer(group1Poly2, dist = 0)
      group1Poly3=sf::st_buffer(group1Poly3, dist = 0)
      
      group1Poly = st_union(group1Poly1,group1Poly2)
      group1Poly = st_sf(st_union(group1Poly,group1Poly3))
      
      grp1fin = st_intersection(grp1,group1Poly)
      grp1fin=sf::st_buffer(grp1fin, dist = 0)
      st_geometry(grp1fin) <- "overlop"
      
      group2Poly1=sf::st_buffer(group2Poly1, dist = 0)
      group2Poly2=sf::st_buffer(group2Poly2, dist = 0)
      group2Poly3=sf::st_buffer(group2Poly3, dist = 0)
      
      group2Poly = st_union(group2Poly1,group2Poly2)
      group2Poly = st_sf(st_union(group2Poly,group2Poly3))%>% 
        st_collection_extract(c("POLYGON"))
      
      grp2fin = st_intersection(grp2,group2Poly)
      grp2fin=sf::st_buffer(grp2fin, dist = 0)
      st_geometry(grp2fin) <- "overlop"
      
      group3Poly1=sf::st_buffer(group3Poly1, dist = 0)
      group3Poly2=sf::st_buffer(group3Poly2, dist = 0)
      group3Poly3=sf::st_buffer(group3Poly3, dist = 0)
      
      group3Poly = st_union(group3Poly1,group3Poly2)
      group3Poly = st_sf(st_union(group3Poly,group3Poly3))
      
      grp3fin = st_intersection(grp3,group3Poly)%>% 
        st_collection_extract(c("POLYGON"))
      grp3fin=sf::st_buffer(grp3fin, dist = 0)
      st_geometry(grp3fin) <- "overlop"
      
      
      if('POLYGON' %in% st_geometry_type(st_difference(grp1fin,grp3fin), by_geometry = T)|
         'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp1fin,grp3fin), by_geometry = T)){
        Polyunion1 =st_difference(st_union(grp1fin),st_union(grp3fin))%>% 
          st_collection_extract(c("POLYGON"))
        Polyunion1 =st_intersection(grp1fin,Polyunion1)%>% 
          st_collection_extract(c("POLYGON"))
        st_geometry(Polyunion1) <- "overlop"
        if('POLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)|
           'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)){
          Polyunion2 =st_difference(st_union(grp2fin),st_union(st_union(grp1fin),st_union(grp3fin)))%>% 
            st_collection_extract(c("POLYGON"))
          Polyunion2 =st_intersection(grp2fin,Polyunion2)%>% 
            st_collection_extract(c("POLYGON"))
          st_geometry(Polyunion2) <- "overlop"
        }else{
          Polyunion2=grp2fin
        }
        Polyunion =rbind(Polyunion1,Polyunion2,grp3fin) 
      }else{
        Polyunion1=grp1fin
        if('POLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)|
           'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)){
          Polyunion2 =st_difference(st_union(grp2fin),st_union(st_union(grp1fin),st_union(grp3fin)))%>% 
            st_collection_extract(c("POLYGON"))
          Polyunion2 =st_intersection(grp2fin,Polyunion2)%>% 
            st_collection_extract(c("POLYGON"))
          st_geometry(Polyunion2) <- "overlop"
          Polyunion =rbind(Polyunion1,Polyunion2,grp3fin) 
        }else{
          Polyunion=rbind(grp1fin,grp2fin,grp3fin) 
        }
      }
      
    }
    if(local[4,7]=="LD"){
      
      group1Poly1 <- CreatePolygon(c(RULD_PS[1],RULD_PS[2]), c(RU$x,RU$y), c(RULU_PS[1],RULU_PS[2]),c(RULD_PS[1],RULD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group1Poly2 <- CreatePolygon(c(RDLD_PS[1],RDLD_PS[2]), c(RD$x,RD$y), c(RDLU_PS[1],RDLU_PS[2]),c(RDLD_PS[1],RDLD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group1Poly3 <- CreatePolygon(c(LDRD_PS[1],LDRD_PS[2]), c(LD$x,LD$y), c(LDRU_PS[1],LDRU_PS[2]),c(LURU_PS[1],LURU_PS[2]),c(LU$x,LU$y),c(LURD_PS[1],LURD_PS[2]),c(LDRD_PS[1],LDRD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly1 <- CreatePolygon(c(LURD_PS[1],LURD_PS[2]), c(LU$x,LU$y), c(LULD_PS[1],LULD_PS[2]),c(LURD_PS[1],LURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly2 <- CreatePolygon(c(RURD_PS[1],RURD_PS[2]), c(RU$x,RU$y), c(RULD_PS[1],RULD_PS[2]),c(RURD_PS[1],RURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly3 <- CreatePolygon(c(LDLU_PS[1],LDLU_PS[2]), c(LD$x,LD$y), c(LDRU_PS[1],LDRU_PS[2]),c(RDRU_PS[1],RDRU_PS[2]),c(RD$x,RD$y),c(RDLU_PS[1],RDLU_PS[2]),c(LDLU_PS[1],LDLU_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly1 <- CreatePolygon(c(LULD_PS[1],LULD_PS[2]), c(LU$x,LU$y), c(LURU_PS[1],LURU_PS[2]),c(LULD_PS[1],LULD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly2 <- CreatePolygon(c(RDLD_PS[1],RDLD_PS[2]), c(RD$x,RD$y), c(RDRU_PS[1],RDRU_PS[2]),c(RDLD_PS[1],RDLD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly3 <- CreatePolygon(c(RURD_PS[1],RURD_PS[2]), c(RU$x,RU$y), c(RULU_PS[1],RULU_PS[2]),c(LDLU_PS[1],LDLU_PS[2]),c(LD$x,LD$y),c(LDRD_PS[1],LDRD_PS[2]),c(RURD_PS[1],RURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      
      
      
      group1Poly1=sf::st_buffer(group1Poly1, dist = 0)
      group1Poly2=sf::st_buffer(group1Poly2, dist = 0)
      group1Poly3=sf::st_buffer(group1Poly3, dist = 0)
      
      group1Poly = st_union(group1Poly1,group1Poly2)
      group1Poly = st_sf(st_union(group1Poly,group1Poly3))
      
      grp1fin = st_intersection(grp1,group1Poly)
      grp1fin=sf::st_buffer(grp1fin, dist = 0)
      st_geometry(grp1fin) <- "overlop"
      
      group2Poly1=sf::st_buffer(group2Poly1, dist = 0)
      group2Poly2=sf::st_buffer(group2Poly2, dist = 0)
      group2Poly3=sf::st_buffer(group2Poly3, dist = 0)
      
      group2Poly = st_union(group2Poly1,group2Poly2)
      group2Poly = st_sf(st_union(group2Poly,group2Poly3))%>% 
        st_collection_extract(c("POLYGON"))
      
      grp2fin = st_intersection(grp2,group2Poly)
      grp2fin=sf::st_buffer(grp2fin, dist = 0)
      st_geometry(grp2fin) <- "overlop"
      
      group3Poly1=sf::st_buffer(group3Poly1, dist = 0)
      group3Poly2=sf::st_buffer(group3Poly2, dist = 0)
      group3Poly3=sf::st_buffer(group3Poly3, dist = 0)
      
      group3Poly = st_union(group3Poly1,group3Poly2)
      group3Poly = st_sf(st_union(group3Poly,group3Poly3))
      
      grp3fin = st_intersection(grp3,group3Poly)%>% 
        st_collection_extract(c("POLYGON"))
      grp3fin=sf::st_buffer(grp3fin, dist = 0)
      st_geometry(grp3fin) <- "overlop"
      
      
      if('POLYGON' %in% st_geometry_type(st_difference(grp1fin,grp3fin), by_geometry = T)|
         'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp1fin,grp3fin), by_geometry = T)){
        Polyunion1 =st_difference(st_union(grp1fin),st_union(grp3fin))%>% 
          st_collection_extract(c("POLYGON"))
        Polyunion1 =st_intersection(grp1fin,Polyunion1)%>% 
          st_collection_extract(c("POLYGON"))
        st_geometry(Polyunion1) <- "overlop"
        if('POLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)|
           'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)){
          Polyunion2 =st_difference(st_union(grp2fin),st_union(st_union(grp1fin),st_union(grp3fin)))%>% 
            st_collection_extract(c("POLYGON"))
          Polyunion2 =st_intersection(grp2fin,Polyunion2)%>% 
            st_collection_extract(c("POLYGON"))
          st_geometry(Polyunion2) <- "overlop"
        }else{
          Polyunion2=grp2fin
        }
        Polyunion =rbind(Polyunion1,Polyunion2,grp3fin) 
      }else{
        Polyunion1=grp1fin
        if('POLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)|
           'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)){
          Polyunion2 =st_difference(st_union(grp2fin),st_union(st_union(grp1fin),st_union(grp3fin)))%>% 
            st_collection_extract(c("POLYGON"))
          Polyunion2 =st_intersection(grp2fin,Polyunion2)%>% 
            st_collection_extract(c("POLYGON"))
          st_geometry(Polyunion2) <- "overlop"
          Polyunion =rbind(Polyunion1,Polyunion2,grp3fin) 
        }else{
          Polyunion=rbind(grp1fin,grp2fin,grp3fin) 
        }
      }
      
    }
    if(local[4,7]=="LU"){
      
      group1Poly1 <- CreatePolygon(c(RULD_PS[1],RULD_PS[2]), c(RU$x,RU$y), c(RULU_PS[1],RULU_PS[2]),c(RULD_PS[1],RULD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group1Poly2 <- CreatePolygon(c(RDLD_PS[1],RDLD_PS[2]), c(RD$x,RD$y), c(RDLU_PS[1],RDLU_PS[2]),c(RDLD_PS[1],RDLD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group1Poly3 <- CreatePolygon(c(LURD_PS[1],LURD_PS[2]), c(LU$x,LU$y), c(LURU_PS[1],LURU_PS[2]),c(LDRU_PS[1],LDRU_PS[2]),c(LD$x,LD$y),c(LDRD_PS[1],LDRD_PS[2]),c(LURD_PS[1],LURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly1 <- CreatePolygon(c(LDLU_PS[1],LDLU_PS[2]), c(LD$x,LD$y), c(LDRU_PS[1],LDRU_PS[2]),c(LDLU_PS[1],LDLU_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly2 <- CreatePolygon(c(RDRU_PS[1],RDRU_PS[2]), c(RD$x,RD$y), c(RDLU_PS[1],RDLU_PS[2]),c(RDRU_PS[1],RDRU_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly3 <- CreatePolygon(c(LURD_PS[1],LURD_PS[2]), c(LU$x,LU$y), c(LULD_PS[1],LULD_PS[2]),c(RULD_PS[1],RULD_PS[2]),c(RU$x,RU$y),c(RURD_PS[1],RURD_PS[2]),c(LURD_PS[1],LURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly1 <- CreatePolygon(c(LDRD_PS[1],LDRD_PS[2]), c(LD$x,LD$y), c(LDLU_PS[1],LDLU_PS[2]),c(LDRD_PS[1],LDRD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly2 <- CreatePolygon(c(RURD_PS[1],RURD_PS[2]), c(RU$x,RU$y), c(RULU_PS[1],RULU_PS[2]),c(RURD_PS[1],RURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly3 <- CreatePolygon(c(LULD_PS[1],LULD_PS[2]), c(LU$x,LU$y), c(LURU_PS[1],LURU_PS[2]),c(RDRU_PS[1],RDRU_PS[2]),c(RD$x,RD$y),c(RDLD_PS[1],RDLD_PS[2]),c(LULD_PS[1],LULD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      
      group1Poly1=sf::st_buffer(group1Poly1, dist = 0)
      group1Poly2=sf::st_buffer(group1Poly2, dist = 0)
      group1Poly3=sf::st_buffer(group1Poly3, dist = 0)
      
      group1Poly = st_union(group1Poly1,group1Poly2)
      group1Poly = st_sf(st_union(group1Poly,group1Poly3))
      
      grp1fin = st_intersection(grp1,group1Poly)
      grp1fin=sf::st_buffer(grp1fin, dist = 0)
      st_geometry(grp1fin) <- "overlop"
      
      group2Poly1=sf::st_buffer(group2Poly1, dist = 0)
      group2Poly2=sf::st_buffer(group2Poly2, dist = 0)
      group2Poly3=sf::st_buffer(group2Poly3, dist = 0)
      
      group2Poly = st_union(group2Poly1,group2Poly2)
      group2Poly = st_sf(st_union(group2Poly,group2Poly3))%>% 
        st_collection_extract(c("POLYGON"))
      
      grp2fin = st_intersection(grp2,group2Poly)
      grp2fin=sf::st_buffer(grp2fin, dist = 0)
      st_geometry(grp2fin) <- "overlop"
      
      group3Poly1=sf::st_buffer(group3Poly1, dist = 0)
      group3Poly2=sf::st_buffer(group3Poly2, dist = 0)
      group3Poly3=sf::st_buffer(group3Poly3, dist = 0)
      
      group3Poly = st_union(group3Poly1,group3Poly2)
      group3Poly = st_sf(st_union(group3Poly,group3Poly3))
      
      grp3fin = st_intersection(grp3,group3Poly)%>% 
        st_collection_extract(c("POLYGON"))
      grp3fin=sf::st_buffer(grp3fin, dist = 0)
      st_geometry(grp3fin) <- "overlop"
      
      
      if('POLYGON' %in% st_geometry_type(st_difference(grp1fin,grp3fin), by_geometry = T)|
         'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp1fin,grp3fin), by_geometry = T)){
        Polyunion1 =st_difference(st_union(grp1fin),st_union(grp3fin))%>% 
          st_collection_extract(c("POLYGON"))
        Polyunion1 =st_intersection(grp1fin,Polyunion1)%>% 
          st_collection_extract(c("POLYGON"))
        st_geometry(Polyunion1) <- "overlop"
        if('POLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)|
           'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)){
          Polyunion2 =st_difference(st_union(grp2fin),st_union(st_union(grp1fin),st_union(grp3fin)))%>% 
            st_collection_extract(c("POLYGON"))
          Polyunion2 =st_intersection(grp2fin,Polyunion2)%>% 
            st_collection_extract(c("POLYGON"))
          st_geometry(Polyunion2) <- "overlop"
        }else{
          Polyunion2=grp2fin
        }
        Polyunion =rbind(Polyunion1,Polyunion2,grp3fin) 
      }else{
        Polyunion1=grp1fin
        if('POLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)|
           'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)){
          Polyunion2 =st_difference(st_union(grp2fin),st_union(st_union(grp1fin),st_union(grp3fin)))%>% 
            st_collection_extract(c("POLYGON"))
          Polyunion2 =st_intersection(grp2fin,Polyunion2)%>% 
            st_collection_extract(c("POLYGON"))
          st_geometry(Polyunion2) <- "overlop"
          Polyunion =rbind(Polyunion1,Polyunion2,grp3fin) 
        }else{
          Polyunion=rbind(grp1fin,grp2fin,grp3fin) 
        }
      }
      
    }
    if(local[4,7]=="RD"){
      
      group1Poly1 <- CreatePolygon(c(LDRU_PS[1],LDRU_PS[2]), c(LD$x,LD$y), c(LDRD_PS[1],LDRD_PS[2]),c(LDRU_PS[1],LDRU_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group1Poly2 <- CreatePolygon(c(LURD_PS[1],LURD_PS[2]), c(LU$x,LU$y), c(LURU_PS[1],LURU_PS[2]),c(LURD_PS[1],LURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group1Poly3 <- CreatePolygon(c(RULD_PS[1],RULD_PS[2]), c(RU$x,RU$y), c(RULU_PS[1],RULU_PS[2]),c(RDLU_PS[1],RDLU_PS[2]),c(RD$x,RD$y),c(RDLD_PS[1],RDLD_PS[2]),c(RULD_PS[1],RULD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly1 <- CreatePolygon(c(LURD_PS[1],LURD_PS[2]), c(LU$x,LU$y), c(LULD_PS[1],LULD_PS[2]),c(LURD_PS[1],LURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly2 <- CreatePolygon(c(RURD_PS[1],RURD_PS[2]), c(RU$x,RU$y), c(RULD_PS[1],RULD_PS[2]),c(RURD_PS[1],RURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group2Poly3 <- CreatePolygon(c(RDLU_PS[1],RDLU_PS[2]), c(RD$x,RD$y), c(RDRU_PS[1],RDRU_PS[2]),c(LDRU_PS[1],LDRU_PS[2]),c(LD$x,LD$y),c(LDLU_PS[1],LDLU_PS[2]),c(RDLU_PS[1],RDLU_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly2 <- CreatePolygon(c(LDRD_PS[1],LDRD_PS[2]), c(LD$x,LD$y), c(LDLU_PS[1],LDLU_PS[2]),c(LDRD_PS[1],LDRD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly1 <- CreatePolygon(c(RURD_PS[1],RURD_PS[2]), c(RU$x,RU$y), c(RULU_PS[1],RULU_PS[2]),c(RURD_PS[1],RURD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly3 <- CreatePolygon(c(LULD_PS[1],LULD_PS[2]), c(LU$x,LU$y), c(LURU_PS[1],LURU_PS[2]),c(RDRU_PS[1],RDRU_PS[2]),c(RD$x,RD$y),c(RDLD_PS[1],RDLD_PS[2]),c(LULD_PS[1],LULD_PS[2]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      
      
      
      group1Poly1=sf::st_buffer(group1Poly1, dist = 0)
      group1Poly2=sf::st_buffer(group1Poly2, dist = 0)
      group1Poly3=sf::st_buffer(group1Poly3, dist = 0)
      
      group1Poly = st_union(group1Poly1,group1Poly2)
      group1Poly = st_sf(st_union(group1Poly,group1Poly3))
      
      grp1fin = st_intersection(grp1,group1Poly)
      grp1fin=sf::st_buffer(grp1fin, dist = 0)
      st_geometry(grp1fin) <- "overlop"
      
      group2Poly1=sf::st_buffer(group2Poly1, dist = 0)
      group2Poly2=sf::st_buffer(group2Poly2, dist = 0)
      group2Poly3=sf::st_buffer(group2Poly3, dist = 0)
      
      group2Poly = st_union(group2Poly1,group2Poly2)
      group2Poly = st_sf(st_union(group2Poly,group2Poly3))%>% 
        st_collection_extract(c("POLYGON"))
      
      grp2fin = st_intersection(grp2,group2Poly)
      grp2fin=sf::st_buffer(grp2fin, dist = 0)
      st_geometry(grp2fin) <- "overlop"
      
      group3Poly1=sf::st_buffer(group3Poly1, dist = 0)
      group3Poly2=sf::st_buffer(group3Poly2, dist = 0)
      group3Poly3=sf::st_buffer(group3Poly3, dist = 0)
      
      group3Poly = st_union(group3Poly1,group3Poly2)
      group3Poly = st_sf(st_union(group3Poly,group3Poly3))
      
      grp3fin = st_intersection(grp3,group3Poly)%>% 
        st_collection_extract(c("POLYGON"))
      grp3fin=sf::st_buffer(grp3fin, dist = 0)
      st_geometry(grp3fin) <- "overlop"
      
      if('POLYGON' %in% st_geometry_type(st_difference(grp1fin,grp3fin), by_geometry = T)|
         'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp1fin,grp3fin), by_geometry = T)){
        Polyunion1 =st_difference(st_union(grp1fin),st_union(grp3fin))%>% 
          st_collection_extract(c("POLYGON"))
        Polyunion1 =st_intersection(grp1fin,Polyunion1)%>% 
          st_collection_extract(c("POLYGON"))
        st_geometry(Polyunion1) <- "overlop"
        if('POLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)|
           'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)){
          Polyunion2 =st_difference(st_union(grp2fin),st_union(st_union(grp1fin),st_union(grp3fin)))%>% 
            st_collection_extract(c("POLYGON"))
          Polyunion2 =st_intersection(grp2fin,Polyunion2)%>% 
            st_collection_extract(c("POLYGON"))
          st_geometry(Polyunion2) <- "overlop"
        }else{
          Polyunion2=grp2fin
        }
        Polyunion =rbind(Polyunion1,Polyunion2,grp3fin) 
      }else{
        Polyunion1=grp1fin
        if('POLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)|
           'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)){
          Polyunion2 =st_difference(st_union(grp2fin),st_union(st_union(grp1fin),st_union(grp3fin)))%>% 
            st_collection_extract(c("POLYGON"))
          Polyunion2 =st_intersection(grp2fin,Polyunion2)%>% 
            st_collection_extract(c("POLYGON"))
          st_geometry(Polyunion2) <- "overlop"
          Polyunion =rbind(Polyunion1,Polyunion2,grp3fin) 
        }else{
          Polyunion=rbind(grp1fin,grp2fin,grp3fin) 
        }
      }
      
    }
  }else{
    
    grp3lin1=ptincerts3V(RURD,LULD)
    grp3lin2=ptincerts3H(LURU,LDRD)
    
    group1Poly1 <- CreatePolygon(c(LU$x,LU$y), c(LURD_PS[1],LURD_PS[2]), c(LURU_PS[1],LURU_PS[2]),c(LU$x,LU$y))%>% 
      list %>% 
      st_polygon %>% 
      st_sfc
    group1Poly2 <- CreatePolygon(c(LU$x,LU$y), c(RDLU_PS[1],RDLU_PS[2]), c(RULU_PS[1],RULU_PS[2]),c(LU$x,LU$y))%>% 
      list %>% 
      st_polygon %>% 
      st_sfc
    group1Poly3 <- CreatePolygon(c(LD$x,LD$y), c(LDRD_PS[1],LDRD_PS[2]), c(LDRU_PS[1],LDRU_PS[2]),c(LD$x,LD$y))%>% 
      list %>% 
      st_polygon %>% 
      st_sfc
    group1Poly4 <- CreatePolygon(c(LD$x,LD$y), c(RDLD_PS[1],RDLD_PS[2]), c(RULD_PS[1],RULD_PS[2]),c(LD$x,LD$y))%>% 
      list %>% 
      st_polygon %>% 
      st_sfc
    
    group2Poly1 <- CreatePolygon(c(RD$x,RD$y), c(RDLU_PS[1],RDLU_PS[2]), c(RDRU_PS[1],RDRU_PS[2]),c(RD$x,RD$y))%>% 
      list %>% 
      st_polygon %>% 
      st_sfc
    group2Poly2 <- CreatePolygon(c(RD$x,RD$y), c(LURD_PS[1],LURD_PS[2]), c(RURD_PS[1],RURD_PS[2]),c(RD$x,RD$y))%>% 
      list %>% 
      st_polygon %>% 
      st_sfc
    group2Poly3 <- CreatePolygon(c(LD$x,LD$y), c(LDLU_PS[1],LDLU_PS[2]), c(LDRU_PS[1],LDRU_PS[2]),c(LD$x,LD$y))%>% 
      list %>% 
      st_polygon %>% 
      st_sfc
    group2Poly4 <- CreatePolygon(c(LD$x,LD$y), c(LULD_PS[1],LULD_PS[2]), c(RULD_PS[1],RULD_PS[2]),c(LD$x,LD$y))%>% 
      list %>% 
      st_polygon %>% 
      st_sfc
    
    group3Poly5 <- CreatePolygon(c(RD$x,RD$y), c(RU$x,RU$y), c(LU$x,LU$y),c(LD$x,LD$y),c(RD$x,RD$y))%>% 
      list %>% 
      st_polygon %>% 
      st_sfc
    
    if(length(grp3lin1)==10){
      group3Poly1 <- CreatePolygon(c(grp3lin1[1],grp3lin1[6]), c(grp3lin1[3],grp3lin1[8]), c(grp3lin1[5],grp3lin1[10]),c(grp3lin1[1],grp3lin1[6]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly2 <- CreatePolygon(c(grp3lin1[1],grp3lin1[6]), c(grp3lin1[2],grp3lin1[7]), c(grp3lin1[4],grp3lin1[9]),c(grp3lin1[1],grp3lin1[6]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
    }else{
      group3Poly1 <- CreatePolygon(c(grp3lin1[2],grp3lin1[6]), c(grp3lin1[3],grp3lin1[7]), c(RU$x,RU$y),c(LU$x,LU$y),c(grp3lin1[2],grp3lin1[6]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly2 <- CreatePolygon(c(grp3lin1[1],grp3lin1[5]), c(grp3lin1[4],grp3lin1[8]), c(LD$x,LD$y),c(RD$x,RD$y),c(grp3lin1[1],grp3lin1[5]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
    }
    if(length(grp3lin2)==10){
      group3Poly3 <- CreatePolygon(c(grp3lin2[1],grp3lin2[6]), c(grp3lin2[3],grp3lin2[8]), c(grp3lin2[5],grp3lin2[10]),c(grp3lin2[1],grp3lin2[6]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly4 <- CreatePolygon(c(grp3lin2[1],grp3lin2[6]), c(grp3lin2[2],grp3lin2[7]), c(grp3lin2[4],grp3lin2[9]),c(grp3lin2[1],grp3lin2[6]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
    }else{
      group3Poly3 <- CreatePolygon(c(grp3lin2[2],grp3lin2[6]), c(grp3lin2[4],grp3lin2[8]), c(LD$x,LD$y),c(LU$x,LU$y),c(grp3lin2[2],grp3lin2[6]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
      group3Poly4 <- CreatePolygon(c(grp3lin2[1],grp3lin2[5]), c(grp3lin2[3],grp3lin2[7]), c(RD$x,RD$y),c(RU$x,RU$y),c(grp3lin2[1],grp3lin2[5]))%>% 
        list %>% 
        st_polygon %>% 
        st_sfc
    }
    
    
    group1Poly1=sf::st_buffer(group1Poly1, dist = 0)
    group1Poly2=sf::st_buffer(group1Poly2, dist = 0)
    group1Poly3=sf::st_buffer(group1Poly3, dist = 0)
    group1Poly4=sf::st_buffer(group1Poly4, dist = 0)
    
    group1Poly = st_union(group1Poly1,group1Poly2)
    group1Poly = st_union(group1Poly,group1Poly3)
    group1Poly = st_sf(st_union(group1Poly,group1Poly4))
    
    grp1fin = st_intersection(grp1,group1Poly)
    grp1fin=sf::st_buffer(grp1fin, dist = 0)
    
    group2Poly1=sf::st_buffer(group2Poly1, dist = 0)
    group2Poly2=sf::st_buffer(group2Poly2, dist = 0)
    group2Poly3=sf::st_buffer(group2Poly3, dist = 0)
    group2Poly4=sf::st_buffer(group2Poly4, dist = 0)
    
    group2Poly = st_union(group2Poly1,group2Poly2)
    group2Poly = st_union(group2Poly,group2Poly3)
    group2Poly = st_sf(st_union(group2Poly,group2Poly4))
    
    grp2fin = st_intersection(grp2,group2Poly)
    grp2fin=sf::st_buffer(grp2fin, dist = 0)
    group3Poly1=sf::st_buffer(group3Poly1, dist = 0)
    group3Poly2=sf::st_buffer(group3Poly2, dist = 0)
    group3Poly3=sf::st_buffer(group3Poly3, dist = 0)
    group3Poly4=sf::st_buffer(group3Poly4, dist = 0)
    group3Poly5=sf::st_buffer(group3Poly5, dist = 0)
    
    group3Poly = st_union(group3Poly1,group3Poly2)
    group3Poly = st_union(group3Poly,group3Poly3)
    group3Poly = st_union(group3Poly,group3Poly4)
    group3Poly = st_sf(st_union(group3Poly,group3Poly5))
    
    grp3fin = st_intersection(grp3,group3Poly)%>% 
      st_collection_extract(c("POLYGON"))
    grp3fin=sf::st_buffer(grp3fin, dist = 0)
    
    st_geometry(grp3fin) <- "overlop"
    
    if('POLYGON' %in% st_geometry_type(st_difference(grp1fin,grp3fin), by_geometry = T)|
       'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp1fin,grp3fin), by_geometry = T)){
      Polyunion1 =st_difference(st_union(grp1fin),st_union(grp3fin))%>% 
        st_collection_extract(c("POLYGON"))
      Polyunion1 =st_intersection(grp1fin,Polyunion1)%>% 
        st_collection_extract(c("POLYGON"))
      st_geometry(Polyunion1) <- "overlop"
      if('POLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)|
         'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)){
        Polyunion2 =st_difference(st_union(grp2fin),st_union(st_union(grp1fin),st_union(grp3fin)))%>% 
          st_collection_extract(c("POLYGON"))
        Polyunion2 =st_intersection(grp2fin,Polyunion2)%>% 
          st_collection_extract(c("POLYGON"))
        st_geometry(Polyunion2) <- "overlop"
      }else{
        Polyunion2=grp2fin
      }
      Polyunion =rbind(Polyunion1,Polyunion2,grp3fin) 
    }else{
      Polyunion1=grp1fin
      if('POLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)|
         'MULTIPOLYGON' %in% st_geometry_type(st_difference(grp2fin,grp3fin), by_geometry = T)){
        Polyunion2 =st_difference(st_union(grp2fin),st_union(st_union(grp1fin),st_union(grp3fin)))%>% 
          st_collection_extract(c("POLYGON"))
        Polyunion2 =st_intersection(grp2fin,Polyunion2)%>% 
          st_collection_extract(c("POLYGON"))
        st_geometry(Polyunion2) <- "overlop"
        Polyunion =rbind(Polyunion1,Polyunion2,grp3fin) 
      }else{
        Polyunion=rbind(grp1fin,grp2fin,grp3fin) 
      }
    }
    
  }
  
  
  Polyunionfin = st_intersection(Polyunion,zoi[i,])%>% 
    st_collection_extract(c("POLYGON"))
  Polyunionfin =unique(Polyunionfin)
  if(i==1){
    ovlap = Polyunionfin
  }else{
    ovlap = rbind(ovlap,Polyunionfin)
  }
  
}
timeend<-Sys.time()
runningtime<-timeend-timestart
print(runningtime)


#####Exporting vector results of point-based UAI
ovlapwt = st_sf(ovlap)
ovlapwt =ovlapwt[,-2]
sf::st_write(ovlapwt, "ovlapwt70.shp")








