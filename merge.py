import pandas as pd

# 输入文件路径
file1 = "path/to/input_file.xlsx"
file2 = "path/to/input_file.xlsx"
output_file = "merged.xlsx"

# 读取两个 Excel 文件的第一个表
df1 = pd.read_excel(file1, sheet_name=0)
df2 = pd.read_excel(file2, sheet_name=0)

# 假设两边要匹配的都是第一列
col1 = df1.columns[0]
col2 = df2.columns[0]

# 按照第一列匹配，把 df2 的全部列带上
merged = pd.merge(df1[[col1]], df2, left_on=col1, right_on=col2, how="left")

# 保存结果到新 Excel
merged.to_excel(output_file, index=False)

print("合并完成，输出文件为:", output_file)
