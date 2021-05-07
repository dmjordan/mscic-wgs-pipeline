import hail as hl
import sysconfig

hl.init()

print("Successfully initialized Hail")
print("Spark configuration:")
for key, value in hl.spark_context().getConf().getAll():
    print(f"\t{key}: {value}")
print("Python configuration:")
for key, value in sysconfig.get_config_vars().items():
    print(f"\t{key}: {value}")

