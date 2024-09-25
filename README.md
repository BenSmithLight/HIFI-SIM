# 主要函数

主要可执行函数包括4类：

1. `HIFI_from_image_stack_main.m`：从包含多帧的单个tif中执行HIFI-SIM
2. `HIFI_from_Multiple_single_images_main.m`：从多个单帧tif中执行HIFI-SIM
3. `HIFI_from_MRC`：从MRC文件中执行HIFI-SIM（缺少主函数，需要自行调用）
4. `HIFI_Main.m`：原版UI程序

---

# 文件夹结构

注意：

1. `rawdata`文件夹与`HIFI-SIM`文件夹同级

- HIFI-SIM
    - HIFI_from_image_stack_main.m
    - HIFI_from_Multiple_single_images_main.m
    - ……
- rawdata
    - image1
        - xxx.tif
- processed_data
    - image1
        - xxx_HIFI
            - output.tif
        - xxx_merge
            - merge.tif