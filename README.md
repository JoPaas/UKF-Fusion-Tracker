## Unscented Kalman Filter (UKF) Fusion Tracker
Johannes Paas

To execute the code to the following from root level:
1. mkdir build
2. cd build
3. cmake ..
4. make
5. ./UnscentedKF

[//]: # (Image References)
[image1]: ./initial_guess_3_003.png
[image2]: ./educated_guess_2_08.png
[image3]: ./educated_guess_2_08_D2.png

---

### Initial attempt
```python
std_a: 3.0 m/s²
std_yawdd: 0.03 rad/s²

RMSE
0.3530
0.3289
0.7803
0.7390
```
![alt_text][image1]

---

### Final Result Dataset_1
```python
std_a: 2.0 m/s²
std_yawdd: 0.8 rad/s²

RMSE
0.0704
0.0820
0.3292
0.2384
```
![alt_text][image2]

The initial guess for process noise was too high, since a bycicle does not accelerate or turn that hard. After dropping the Values in a couple of experiments I settled with 2 m/s² for the acceleration standard deviation, which means that the bicycle accelerates with less than 4 m/s² in 95% of the cases and 0.8 rad/s² for the yaw accelerations standard deviateion.

---

### Final Result Dataset_2
```python
std_a: 2.0 m/s²
std_yawdd: 0.8 rad/s²

RMSE
0.0705
0.0696
0.1788
0.2202
```
![alt_text][image3]

Note: The RMSE is smaller for Dataset_2. This might be because the standard deviation is larger at the beginning and can be further investigated.

---

### Final Thoughts

The best result with the previous EKF-approach was:
```python
RMSE
0.0981
0.0858
0.4659
0.4755
```

Since the simulated data was the same the RMSE proves, that the motion model used in the UKF gives an advantage when tracking a vehicle. The knowledge of the way the object moves improves especially the RMSE of vx and vy (by ~50%), which is very important for behavior prediction in a autonomous vehicle.

