Please set the working directory to the folder of `All_lipids_plots-v2.R` in RStudio or use `setwd` 

To save the image using `All_lipids_plots-v2.R`:
- use RStudio and run the code by segments, then save image using RStudio
- add `ggsave('./img/img_name.png')` after `ggplot` function .e.g.
    - ```
      ggplot(data = lpc_df, mapping = aes(x = RT_P, y = `kmd(h)`, color = Bulk_db))  + geom_point(aes(shape=Bulk_l), size = 4) +
      scale_shape_manual(values = c(19,17, 15))+
      scale_color_manual(values = cols, aesthetics = c("colour", "fill"),  breaks= c(0, 1,2,3,4,5,6,7,8,9,10,11,12,13)) +
      scale_y_continuous(sec.axis = sec_axis(trans=trans, name="# number of carbons", breaks = c(14, 15, 16, 17, 18, 19, 20, 22, 24))) +
      labs(title = "LPC RT information", x = "Retention Time", y = "KMD(H)") +
      theme_bw()
      # add following line to save image
      ggsave('./img/img_name.png')
      ```
    