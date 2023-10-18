#' Title swap noise function 1
#' @description The swap noise function1, swap noise, in the ZINB or ZIP
#' @param lambda The probability of the sample() function process.
#' @param spots The number of spots.
#' @param se The number of  spatially variable genes (SVGs).
#' @param ns  The number of  non-SVGs.
#' @param type A character.'ZINB' or 'ZIP. The default value is ZINB.
#' @param se.p A number. The zero proportion of zero generation process of
#' SVGs in the streak area.
#' @param se.size  The size para of SVGs in the NB distribution.
#' @param se.mu  For SVGs, the lambda para in the poisson distribution
#' or the mu para in the NB distribution.
#' @param ns.p  A number. The zero proportion of zero generation process of
#' non-SVGs and SVGs in the non-streak area.
#' @param ns.size For non-SVGs and SVGs in the non-streak area,
#' the size para in the NB distribution.
#' @param ns.mu For non-SVGs and SVGs in the non-streak area, the lambda para in the poisson distribution
#' or the mu para in the NB distribution.
#' @param ptn The file name of the pattern png.
#' @param png_dir The dir of the png files.
#' @param noise swap porpotion.
#' @return A data.frame of gene expression counts and coordinates.
#' @export
noise.zi=function(lambda = 0.7, spots, se, ns, type, se.p, se.size = 0.2,
                  se.mu, ns.p, ns.size = 0.2, ns.mu, ptn,png_dir,
                  noise=0.1) {
  ## noise是比例
  win = ceiling(sqrt(spots/lambda))
  coor = spatstat.random::rpoispp(lambda = lambda, win = spatstat.geom::owin(c(0,win), c(0, win)))
  coor.dt = data.frame(num = paste0("c-", 1:length(coor$x)),
                       row = round(coor$x)+1, col = round(coor$y)+1,
                       row.names = paste0("c_",1:length(coor$x)))

  ## 对于灰度图像的处理
  image <- imager::load.image(paste0(png_dir,ptn))  ## png_dir要提供
  # gray_img=imager::grayscale(image)
  re_img=imager::imresize(imager::grayscale(image) ,scale =round((win)/700,2)) ## 先灰度再缩放
  img_coor=to_dt(round(as.matrix(re_img),1))
  img_coor=img_coor[img_coor$value!=0.5,]; img_coor$value[img_coor$value<0.5]<-0;img_coor$value[img_coor$value>0.5]<-1

  ## 合并，并区分marked area
  coor.s1 = merge(coor.dt,img_coor,by=c('row','col'))
  # coor.s1=coor.s1[coor.s1$value!=0.5,]; coor.s1$value[coor.s1$value<0.5]<-0;coor.s1$value[coor.s1$value>0.5]<-1

  ##
  coor.mark=coor.s1 %>% dplyr::filter(value==1) %>% dplyr::select(num,row,col)
  coor.random=coor.s1 %>% dplyr::filter(value==0) %>% dplyr::select(num,row,col)

  exp.mark = matrix(NA, nrow = nrow(coor.mark), ncol = se,
                    dimnames = list(NULL, paste0("se.", 1:se)))
  exp.mark = apply(exp.mark, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.mark),
                zi.p = se.p, size = se.size, mu = se.mu)
    z
  })
  exp.random = matrix(NA, nrow = nrow(coor.random), ncol = se,
                      dimnames = list(NULL, paste0("se.", 1:se)))
  exp.random = apply(exp.random, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.random),
                zi.p = ns.p, size = ns.size, mu = ns.mu)
    z
  })

  ## 交换部分
  swap_od1=sample(1:nrow(exp.mark),size=round(noise*nrow(exp.mark)))
  swap_od2=sample(1:nrow(exp.random),size=round(noise*nrow(exp.mark)))

  swap1=exp.mark[swap_od1,]
  swap2=exp.random[swap_od2,]
  exp.mark[swap_od1,]<-swap2
  exp.random[swap_od2,]<-swap1
  ## END 交换部分



  exp.svg = rbind(exp.mark, exp.random)
  non.coor = rbind(coor.mark, coor.random)
  exp.non = matrix(NA, nrow = nrow(non.coor), ncol = ns,
                   dimnames = list(NULL,paste0("ns.", 1:ns)))
  exp.non = apply(exp.non, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.non),
                zi.p = ns.p, size = ns.size, mu = ns.mu)
    z
  })
  all = cbind(non.coor[c("row", "col")], exp.svg, exp.non)
  # end = all
  # end
  all
}




## 对数正态产生noise,差值很大 ------------------------------

#' Title swap noise function 2
#' @description The swap noise function2, swap noise, in the  logNormal with cut.
#' @param cut cur ,the threshold.
#' @param lambda The probability of the sample() function process.
#' @param spots The number of spots.
#' @param se The number of  spatially variable genes (SVGs).
#' @param ns  The number of  non-SVGs.
#' @param type A character.'ZINB' or 'ZIP. The default value is ZINB.
#' @param se.p A number. The zero proportion of zero generation process of
#' SVGs in the streak area.
#' @param se.size  The size para of SVGs in the NB distribution.
#' @param se.mu  For SVGs, the lambda para in the poisson distribution
#' or the mu para in the NB distribution.
#' @param ns.p  A number. The zero proportion of zero generation process of
#' non-SVGs and SVGs in the non-streak area.
#' @param ns.size For non-SVGs and SVGs in the non-streak area,
#' the size para in the NB distribution.
#' @param ns.mu For non-SVGs and SVGs in the non-streak area, the lambda para in the poisson distribution
#' or the mu para in the NB distribution.
#' @param ptn The file name of the pattern png.
#' @param png_dir The dir of the png files.
#' @param noise swap porpotion.
#' @return A data.frame of gene expression counts and coordinates.
#' @export
noise.norm=function (noise=0.2,cut=5,lambda = 0.7, spots, se, ns, se.p, ns.p, se.sd = 0.2,
                     se.mu, ns.sd = 0.2, ns.mu, ptn, png_dir)
{
  win = ceiling(sqrt(spots/lambda))
  coor = spatstat.random::rpoispp(lambda = lambda, win = spatstat.geom::owin(c(0,
                                                                               win), c(0, win)))
  coor.dt = data.frame(num = paste0("c-", 1:length(coor$x)),
                       row = round(coor$x) + 1, col = round(coor$y) + 1, row.names = paste0("c_",
                                                                                            1:length(coor$x)))
  image <- imager::load.image(paste0(png_dir, ptn))
  re_img = imager::imresize(imager::grayscale(image), scale = round((win)/700,
                                                                    2))
  img_coor = to_dt(round(as.matrix(re_img), 1))
  img_coor = img_coor[img_coor$value != 0.5, ]
  img_coor$value[img_coor$value < 0.5] <- 0
  img_coor$value[img_coor$value > 0.5] <- 1
  coor.s1 = merge(coor.dt, img_coor, by = c("row", "col"))
  coor.mark = coor.s1 %>% dplyr::filter(value == 1) %>% dplyr::select(num,
                                                                      row, col)
  coor.random = coor.s1 %>% dplyr::filter(value == 0) %>% dplyr::select(num,
                                                                        row, col)
  exp.mark = matrix(NA, nrow = nrow(coor.mark), ncol = se,
                    dimnames = list(NULL, paste0("se.", 1:se)))


  exp.mark = apply(exp.mark, 2, function(y) {
    z = simu_part_norm(subject.n = nrow(exp.mark), zi.p = se.p, cut=cut,
                       sd = se.sd, mu = se.mu)
    z
  })
  exp.random = matrix(NA, nrow = nrow(coor.random), ncol = se,
                      dimnames = list(NULL, paste0("se.", 1:se)))
  exp.random = apply(exp.random, 2, function(y) {
    z = simu_norm(subject.n = nrow(exp.random), zi.p = ns.p,
                  sd = ns.sd, mu = ns.mu)
    z
  })

  ## noise !! 交换部分
  swap_od1=sample(1:nrow(exp.mark),size=round(noise*nrow(exp.mark)))
  swap_od2=sample(1:nrow(exp.random),size=round(noise*nrow(exp.mark)))

  swap1=exp.mark[swap_od1,]
  swap2=exp.random[swap_od2,]
  exp.mark[swap_od1,]<-swap2
  exp.random[swap_od2,]<-swap1
  ## END 交换部分

  exp.svg = rbind(exp.mark, exp.random)
  non.coor = rbind(coor.mark, coor.random)
  exp.non = matrix(NA, nrow = nrow(non.coor), ncol = ns, dimnames = list(NULL,
                                                                         paste0("ns.", 1:ns)))
  exp.non = apply(exp.non, 2, function(y) {
    z = simu_norm(subject.n = nrow(exp.non), zi.p = ns.p,
                  sd = ns.sd, mu = ns.mu)
    z
  })
  all = cbind(non.coor[c("row", "col")], exp.svg, exp.non)
  end = all
  end
}





## 1004新noise,差值很大 ------------------------------
#' Title swap noise function 3
#' @description The swap noise function3, swap noise, in the ZINB or ZIP.
#' The diff between SVG and non-SVG is large.
#' @param lambda The probability of the sample() function process.
#' @param spots The number of spots.
#' @param se The number of  spatially variable genes (SVGs).
#' @param ns  The number of  non-SVGs.
#' @param type A character.'ZINB' or 'ZIP. The default value is ZINB.
#' @param se.p A number. The zero proportion of zero generation process of
#' SVGs in the streak area.
#' @param se.size  The size para of SVGs in the NB distribution.
#' @param se.mu  For SVGs, the lambda para in the poisson distribution
#' or the mu para in the NB distribution.
#' @param ns.p  A number. The zero proportion of zero generation process of
#' non-SVGs and SVGs in the non-streak area.
#' @param ns.size For non-SVGs and SVGs in the non-streak area,
#' the size para in the NB distribution.
#' @param ns.mu For non-SVGs and SVGs in the non-streak area, the lambda para in the poisson distribution
#' or the mu para in the NB distribution.
#' @param ptn The file name of the pattern png.
#' @param png_dir The dir of the png files.
#' @param noise swap porpotion.
#' @return A data.frame of gene expression counts and coordinates.
#' @export

noise.zi2=function(lambda = 0.7, spots, se, ns, type, se.p, se.size = 0.2,
                   se.mu, ns.p, ns.size = 0.2, ns.mu, ptn,png_dir,
                   noise=0.1) {
  ## noise是比例
  win = ceiling(sqrt(spots/lambda))
  coor = spatstat.random::rpoispp(lambda = lambda, win = spatstat.geom::owin(c(0,win), c(0, win)))
  coor.dt = data.frame(num = paste0("c-", 1:length(coor$x)),
                       row = round(coor$x)+1, col = round(coor$y)+1,
                       row.names = paste0("c_",1:length(coor$x)))

  ## 对于灰度图像的处理
  image <- imager::load.image(paste0(png_dir,ptn))  ## png_dir要提供
  # gray_img=imager::grayscale(image)
  re_img=imager::imresize(imager::grayscale(image) ,scale =round((win)/700,2)) ## 先灰度再缩放
  img_coor=to_dt(round(as.matrix(re_img),1))
  img_coor=img_coor[img_coor$value!=0.5,]; img_coor$value[img_coor$value<0.5]<-0;img_coor$value[img_coor$value>0.5]<-1

  ## 合并，并区分marked area
  coor.s1 = merge(coor.dt,img_coor,by=c('row','col'))
  # coor.s1=coor.s1[coor.s1$value!=0.5,]; coor.s1$value[coor.s1$value<0.5]<-0;coor.s1$value[coor.s1$value>0.5]<-1

  ##
  coor.mark=coor.s1 %>% dplyr::filter(value==1) %>% dplyr::select(num,row,col)
  coor.random=coor.s1 %>% dplyr::filter(value==0) %>% dplyr::select(num,row,col)

  exp.mark = matrix(NA, nrow = nrow(coor.mark), ncol = se,
                    dimnames = list(NULL, paste0("se.", 1:se)))
  exp.mark = apply(exp.mark, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.mark),
                zi.p = se.p, size = se.size, mu = se.mu)
    z
  })
  exp.mark=exp.mark+4    # 保证远高于 non-merker area
  exp.random = matrix(NA, nrow = nrow(coor.random), ncol = se,
                      dimnames = list(NULL, paste0("se.", 1:se)))
  exp.random = apply(exp.random, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.random),
                zi.p = ns.p, size = ns.size, mu = ns.mu)
    z
  })

  ## noise !! 交换部分
  swap_od1=sample(1:nrow(exp.mark),size=round(noise*nrow(exp.mark)))
  swap_od2=sample(1:nrow(exp.random),size=round(noise*nrow(exp.mark)))

  swap1=exp.mark[swap_od1,]
  swap2=exp.random[swap_od2,]
  exp.mark[swap_od1,]<-swap2
  exp.random[swap_od2,]<-swap1
  ## END 交换部分



  exp.svg = rbind(exp.mark, exp.random)
  non.coor = rbind(coor.mark, coor.random)
  exp.non = matrix(NA, nrow = nrow(non.coor), ncol = ns,
                   dimnames = list(NULL,paste0("ns.", 1:ns)))
  exp.non = apply(exp.non, 2, function(y) {
    z = simu_zi(family = type, subject.n = nrow(exp.non),
                zi.p = ns.p, size = ns.size, mu = ns.mu)
    z
  })
  all = cbind(non.coor[c("row", "col")], exp.svg, exp.non)
  # end = all
  # end
  all
}




