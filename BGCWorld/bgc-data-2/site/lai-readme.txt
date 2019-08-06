先跑-m ./ini/start.ini 得到每天模型自己输出的lai。

将qh_out.dayout.ascii最后一列的输出拷贝到excel，修改需要改动的lai值，其中填充值为-1.制作调整后的lai强制输入文件为lai.txt

再次运行-m ./ini/start.ini  -lai C:\Users\xialang2012\Documents\BGCWorld\BGCWorld\bgc-data-2\site\lai.txt

比较两次的年度输出结果