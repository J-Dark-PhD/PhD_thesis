import numpy as np

max_size_range = []
test_list = []

component_to_test = "everything"

# optimised values 0.5%
optimised_5_value = 5.47e-04
optimised_1_value = 3.23e-04

# optimised ref values
ref_lipb = 4.405790775123394e21
ref_structure = 2.039960510562691e20
ref_baffle = 1.0916057608814874e18
ref_tungsten = 3.5920093644115635e17
ref_pipes = 9.45769694335509e18

optimised_lipb = 7.77e-4

test_list = []
if component_to_test == "everything":
    # define max range
    inital_value = 0.005
    max_size_range.append(inital_value)
    for i in range(1, 28):
        max_size_range.append(inital_value * 0.9**i)
    # for value in max_size_range:
    #     print("{:.2e}".format(value))
    # print(max_size_range)

    # define cell sizes in each volume
    for size in max_size_range:
        lipb_size = size
        structure_size = size
        baffle_plate_size = size
        tungsten_size = size
        bz_pipes_size = size

    test_list.append(
        [
            lipb_size,
            structure_size,
            baffle_plate_size,
            tungsten_size,
            bz_pipes_size,
        ]
    )

elif component_to_test == "lipb":
    # define max cell size range
    no_cases = 26
    max_size_range = []
    for i in range(0, no_cases):
        max_size_range.append(np.round(optimised_1_value * 1.05**i, decimals=6))
    # for value in max_size_range:
    #     print("{:.2e}".format(value))

    # define sizes in each volume
    for size in max_size_range:
        lipb_size = size
        structure_size = optimised_1_value
        baffle_plate_size = optimised_1_value
        tungsten_size = optimised_1_value
        bz_pipes_size = optimised_1_value
        test_list.append(
            [
                lipb_size,
                structure_size,
                baffle_plate_size,
                tungsten_size,
                bz_pipes_size,
            ]
        )

elif component_to_test == "structure":
    # define max cell size range
    no_cases = 26
    max_size_range = []
    for i in range(0, no_cases):
        max_size_range.append(np.round(optimised_1_value * 1.05**i, decimals=6))
    # for value in max_size_range:
    #     print("{:.2e}".format(value))

    # define sizes in each volume
    for size in max_size_range:
        lipb_size = optimised_lipb
        structure_size = size
        baffle_plate_size = optimised_1_value
        tungsten_size = optimised_1_value
        bz_pipes_size = optimised_1_value
        test_list.append(
            [
                lipb_size,
                structure_size,
                baffle_plate_size,
                tungsten_size,
                bz_pipes_size,
            ]
        )

else:
    print("component not recognised")
    quit()

test_list = np.array(test_list)
max_size_range = np.array(max_size_range)
# print(test_list[0])
