// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Extended Kalman Filter
//
// ####### STUDENT FILE #######
//
// Usage:
// -Rename this file to "kalmanfilter.cpp" if you want to use this code.

#include "kalmanfilter.h"
#include "utils.h"

// -------------------------------------------------- //
// YOU CAN USE AND MODIFY THESE CONSTANTS HERE
constexpr double ACCEL_STD = 1.0;
constexpr double GYRO_STD = 0.01/180.0 * M_PI;
constexpr double INIT_VEL_STD = 10.0;
constexpr double INIT_PSI_STD = 45.0/180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;
constexpr double LIDAR_RANGE_STD = 3.0;
constexpr double LIDAR_THETA_STD = 0.02;
constexpr double GYRO_BIAS_STD = 0.0;

constexpr int num_state = 5;
constexpr int num_z = 2;
// -------------------------------------------------- //

void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map)
{
    // Assume No Correlation between the Measurements and Update Sequentially
    for(const auto& meas : dataset) {handleLidarMeasurement(meas, map);}
}

int matchLandmarks(LidarMeasurement meas, const BeaconMap& map, VectorXd state, int true_id=NULL)
{
    VectorXd z = Vector2d::Zero();
    z << meas.range,meas.theta;

    double x_b = meas.range*cos(meas.theta);
    double y_b = meas.range*sin(meas.theta);

    std::vector<BeaconData> beacons_within_range = map.getBeaconsWithinRange(state(0), state(1),100.0);
    std::vector<double> landmark_dist;

    if(beacons_within_range.size() == 0){
        return -1;
    }
    for(auto beacon: beacons_within_range)
    {
        // std::cout<<"beacon id: "<<beacon.id<<" pos: "<< beacon.x<<" "<<beacon.y<<std::endl;
        double delta_x = beacon.x - state(0);
        double delta_y = beacon.y - state(1);
        double r = sqrt(delta_x*delta_x + delta_y*delta_y);
        double theta = wrapAngle(atan2(delta_y, delta_x) - state(2));

        double beacon_x = r*cos(theta);
        double beacon_y = r*sin(theta);

        double x_diff = beacon_x - x_b;
        double y_diff = beacon_y - y_b;
        double dist = sqrt(x_diff*x_diff+y_diff*y_diff);
        landmark_dist.push_back(dist);
        std::cout<<"beacons within range: "<< beacon.id<<" dist: "<<dist<<std::endl;
    }
    int min_idx = std::min_element(landmark_dist.begin(), landmark_dist.end())-landmark_dist.begin();

    std::cout<<"Truth: "<< true_id<< " calculated: "<< beacons_within_range[min_idx].id<<std::endl;
    if(landmark_dist[min_idx]<300)
    {
        return beacons_within_range[min_idx].id;
    }else{
        return -1;
    }
}

void KalmanFilter::handleLidarMeasurement(LidarMeasurement meas, const BeaconMap& map)
{   
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Update Step for the Lidar Measurements in the 
        // section below.
        // HINT: use the wrapAngle() function on angular values to always keep angle
        // values within correct range, otherwise strange angle effects might be seen.
        // HINT: You can use the constants: LIDAR_RANGE_STD, LIDAR_THETA_STD
        // HINT: The mapped-matched beacon position can be accessed by the variables
        // map_beacon.x and map_beacon.y
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE

        BeaconData map_beacon = map.getBeaconWithId(meas.id); // Match Beacon with built in Data Association Id
        if (meas.id == -1 || map_beacon.id == -1){
            int id = matchLandmarks(meas, map,state);
            map_beacon = map.getBeaconWithId(id);
        }
        // std::cout<<"meas.id: "<<meas.id<<" beacon.id: "<<map_beacon.id<<std::endl;

        if (meas.id != -1 && map_beacon.id != -1)
        {           
            // The map matched beacon positions can be accessed using: map_beacon.x AND map_beacon.y
            VectorXd z = Vector2d();
            z <<meas.range, meas.theta;
            
            double delta_x = map_beacon.x - state(0);
            double delta_y = map_beacon.y - state(1);
            double r = sqrt(delta_x*delta_x + delta_y*delta_y);
            double theta = atan2(delta_y,delta_x) - state(2);
            theta = wrapAngle(theta);

            VectorXd z_hat = Vector2d();
            z_hat << r, theta;
            VectorXd y = z - z_hat;
            y(1) = wrapAngle(y(1));
            double dist = sqrt(y(0)*y(0)+y(1)*y(1));
            MatrixXd R = Matrix2d();
            R << (LIDAR_RANGE_STD*LIDAR_RANGE_STD), 0, 0,(LIDAR_THETA_STD*LIDAR_THETA_STD);
        
            
            MatrixXd H = MatrixXd(num_z, num_state);
            H(0, 0) = -delta_x/r;
            H(0, 1) = -delta_y/r;
            H(0, 2) = 0;
            H(0, 3) = 0;
            H(0, 4) = 0;
            H(1, 0) = delta_y/r/r;
            H(1, 1) = -delta_x/r/r;
            H(1, 2) = -1;
            H(1, 3) = 0;
            H(1, 4) = 0;

            MatrixXd S = H * cov * H.transpose() + R;
            MatrixXd K = cov * H.transpose() * S.inverse();
            state = state + K * y;
            MatrixXd I = MatrixXd::Identity(num_state, num_state);
            cov = ( I - K*H)*cov;

            // Estimate vehicle heading
            double x_diff = map_beacon.x - state(0);
            double y_diff = map_beacon.y - state(1);
            double psi = wrapAngle(atan2(y_diff, x_diff) - meas.theta);
            state(2) = psi;

        }
               

        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    }
    
}

void KalmanFilter::predictionStep(GyroMeasurement gyro, double dt)
{
    if (isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        // Implement The Kalman Filter Prediction Step for the system in the  
        // section below.
        // HINT: Assume the state vector has the form [PX, PY, PSI, V].
        // HINT: Use the Gyroscope measurement as an input into the prediction step.
        // HINT: You can use the constants: ACCEL_STD, GYRO_STD
        // HINT: use the wrapAngle() function on angular values to always keep angle
        // values within correct range, otherwise strange angle effects might be seen.
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE
        VectorXd x = VectorXd(num_state);
        x <<state(3)*cos(state(2)),\
            state(3)*sin(state(2)),\
            gyro.psi_dot - state(4) + 0,\
            0,\
            0;

        state = state + dt*x;
        state(2) = wrapAngle(state(2));
        MatrixXd F = MatrixXd::Identity(num_state, num_state);
        F(0, 2) = -dt*state(3)*sin(state(2));
        F(0, 3) = dt*cos(state(2));
        F(1, 2) = dt*state(3)*cos(state(2));
        F(1, 3) = dt*sin(state(2));
        F(2, 4) = -dt;
        
        MatrixXd Q = MatrixXd::Zero(num_state, num_state);
        Q(2, 2) = dt*dt*GYRO_STD*GYRO_STD+dt*dt*GYRO_BIAS_STD*GYRO_BIAS_STD;
        Q(2, 4) = -dt*GYRO_BIAS_STD*GYRO_BIAS_STD;
        Q(3, 3) = dt*dt*ACCEL_STD*ACCEL_STD;
        Q(4, 2) = -dt*GYRO_BIAS_STD*GYRO_BIAS_STD;
        Q(4, 4) = GYRO_BIAS_STD*GYRO_BIAS_STD;

        cov = F*cov*F.transpose() + Q;


        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    } 
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    // All this code is the same as the LKF as the measurement model is linear
    // so the EKF update state would just produce the same result.
    if(isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        VectorXd z = Vector2d::Zero();
        MatrixXd H = MatrixXd(num_z,  num_state);
        MatrixXd R = Matrix2d::Zero();

        z << meas.x,meas.y;
        H << 1,0,0,0,0,\
             0,1,0,0,0;
        R(0,0) = GPS_POS_STD*GPS_POS_STD;
        R(1,1) = GPS_POS_STD*GPS_POS_STD;

        VectorXd z_hat = H * state;
        VectorXd y = z - z_hat;
        MatrixXd S = H * cov * H.transpose() + R;
        //Check for faulty measurements
        double eta = y.transpose() * S.inverse() * y;
        double Chi_threshold = 5.99;
        // std::cout<<"eta: "<< eta<<"\nthr: "<<Chi_threshold<<std::endl;
        // if( eta < Chi_threshold)
        // {   
            MatrixXd K = cov*H.transpose()*S.inverse();
            state = state + K*y;
            cov = (MatrixXd::Identity(num_state, num_state) - K*H) * cov;
            setState(state);
            setCovariance(cov);
        // }
        
    }
    else
    {
        VectorXd state = VectorXd::Zero(num_state);
        MatrixXd cov = MatrixXd::Zero(num_state, num_state);

        state(0) = meas.x;
        state(1) = meas.y;
        cov(0,0) = GPS_POS_STD*GPS_POS_STD;
        cov(1,1) = GPS_POS_STD*GPS_POS_STD;
        cov(2,2) = INIT_PSI_STD*INIT_PSI_STD;
        cov(3,3) = INIT_VEL_STD*INIT_VEL_STD;
        cov(4,4) = GYRO_BIAS_STD*GYRO_BIAS_STD;

        setState(state);
        setCovariance(cov);
    } 
             
}

Matrix2d KalmanFilter::getVehicleStatePositionCovariance()
{
    Matrix2d pos_cov = Matrix2d::Zero();
    MatrixXd cov = getCovariance();
    if (isInitialised() && cov.size() != 0){pos_cov << cov(0,0), cov(0,1), cov(1,0), cov(1,1);}
    return pos_cov;
}

VehicleState KalmanFilter::getVehicleState()
{
    if (isInitialised())
    {
        VectorXd state = getState(); // STATE VECTOR [X,Y,PSI,V,...]
        return VehicleState(state[0],state[1],state[2],state[3]);
    }
    return VehicleState();
}

void KalmanFilter::predictionStep(double dt){}
