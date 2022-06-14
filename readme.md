# Fast ECI

This repository holds a test program proving the effectiveness of estimating the transformation from an inertial frame to
a fixed frame using simple z-axes rotation from known rotation.

This approach is extremely useful when numerous transformations are needed within close time proximity to each other.

## Usage

Install the environment

```
conda env create -f environment.yml
```

Run the application
```
python main.py
```

Refer to the help text with the `-h` or `--help` options for more options.
