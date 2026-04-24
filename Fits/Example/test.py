import brick
import numpy as np

azr = brick.AZR("./12C+p.azr")
print(azr.config.get_input_values())
print(azr.config.initial_norms)

print(id(azr.config))


params = [2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5, 2.5]
for i in range(1):
    azr.predict(params)[0]
    print("==========================")

print(id(azr.config))
