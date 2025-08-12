###林分模拟抚育
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

extend_point <- function(A, B, distance = 10000) {
  # 输入：
  # A: 点A的坐标，c(x_A, y_A)
  # B: 点B的坐标，c(x_B, y_B)
  # distance: 延长距离，默认为10000米
  # 输出：点C的坐标，c(x_C, y_C)
  
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

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
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

# 函数：从四个点中找到最大面积三角形的三个点和剩余点，若最大三角形面积大于四边形面积，则返回最小面积三角形
# 输入：points，一个4x2矩阵，每行是一个点的x和y坐标
# 输出：列表，包含max_or_min_triangle_indices（面积最大或最小的三角形三个点行号）和remaining_index（剩余点的行号）
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
local=data1[1:4,]
known.pair=matrix(nrow = 2,ncol = 2)
####Define 3 sets,upper left(LU), lower left (LD), upper right (RU), and lower right (RD)
group1=matrix(nrow = 4,ncol = 2)
group2=matrix(nrow = 4,ncol = 2)
group3=matrix(nrow = 4,ncol = 2)
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
  data3$W[i]=w.f(x0,y0,Near4$x,Near4$y) #function 4
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
        concave_point[,3:6] <- hull_points[1,3:6]#避免出现空洞
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
sf::st_write(ovlapwt, "ovlapwt67.shp")







