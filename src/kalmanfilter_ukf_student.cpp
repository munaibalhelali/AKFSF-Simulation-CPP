// ------------------------------------------------------------------------------- //
// Advanced Kalman Filtering and Sensor Fusion Course - Unscented Kalman Filter
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
constexpr double INIT_VEL_STD = 10;
constexpr double INIT_PSI_STD = 45.0/180.0 * M_PI;
constexpr double GPS_POS_STD = 3.0;
constexpr double LIDAR_RANGE_STD = 3.0;
constexpr double LIDAR_THETA_STD = 0.02;
// -------------------------------------------------- //

// ----------------------------------------------------------------------- //
// USEFUL HELPER FUNCTIONS
VectorXd normaliseState(VectorXd state)
{
    state(2) = wrapAngle(state(2));
    return state;
}

VectorXd normaliseLidarMeasurement(VectorXd meas)
{
    meas(1) = wrapAngle(meas(1));
    return meas;
}

std::vector<VectorXd> generateSigmaPoints(VectorXd state, MatrixXd cov)
{
    std::vector<VectorXd> sigmaPoints;
    size_t n = state.size();

    double k = 3.0 - n;
    MatrixXd sqrtCov = cov.llt().matrixL();
    sigmaPoints.push_back(state);

    for(int iState=0; iState<n; ++iState)
    {
        sigmaPoints.push_back(state + sqrt(n+k)*sqrtCov.col(iState));
        sigmaPoints.push_back(state - sqrt(n+k)*sqrtCov.col(iState));
    }
    

    return sigmaPoints;
}

std::vector<double> generateSigmaWeights(unsigned int numStates)
{   
    std::vector<double> weights;
    double k = 3.0 - numStates;
    double W_0 = k/(numStates+k);  
    weights.push_back(W_0);
    double W_i = 1/(2*(numStates+k));
    for(int i=0; i<2*numStates; ++i)
    {
        weights.push_back(W_i);
    }

    return weights;
}

VectorXd vehicleProcessModel(VectorXd aug_state, double psi_dot, double dt)
{
    VectorXd new_state = VectorXd::Zero(4);
    new_state = aug_state.head(4);

    VectorXd u = VectorXd::Zero(4);
    u << aug_state(3)*cos(aug_state(2)),\
        aug_state(3)*sin(aug_state(2)),\
        psi_dot+aug_state(4),\
        aug_state(5);

    new_state = new_state + dt * u;


    return new_state;
}

VectorXd lidarMeasurementModel(VectorXd aug_state, double beaconX, double beaconY)
{
    double delta_x = beaconX - aug_state(0);
    double delta_y = beaconY- aug_state(1);

    VectorXd z_hat = VectorXd::Zero(2);
    z_hat(0) = sqrt(delta_x*delta_x+delta_y*delta_y)+aug_state(4);
    z_hat(1) = wrapAngle(atan2(delta_y, delta_x) - aug_state(2) + aug_state(5));

    return z_hat;
}
// ----------------------------------------------------------------------- //

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
        // Hint: You can use the constants: LIDAR_RANGE_STD, LIDAR_THETA_STD
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE

        BeaconData map_beacon = map.getBeaconWithId(meas.id); // Match Beacon with built in Data Association Id
        if (meas.id != -1 && map_beacon.id != -1) // Check that we have a valid beacon match
        {
            VectorXd z = Vector2d::Zero();
            z << meas.range, meas.theta;

            int num_state = state.size();
            int num_w = 2;
            int num_aug = num_state + num_w;

            MatrixXd R = MatrixXd::Zero(num_w, num_w);
            R(0,0) = LIDAR_RANGE_STD*LIDAR_RANGE_STD;
            R(1,1) = LIDAR_THETA_STD*LIDAR_THETA_STD;

            VectorXd aug_state = VectorXd::Zero(num_aug);
            MatrixXd aug_cov = MatrixXd::Zero(num_aug, num_aug);
            aug_state.head(num_state) = state;
            aug_cov.topLeftCorner(num_state, num_state) = cov;
            aug_cov.bottomRightCorner(num_w, num_w) = R;          

            std::vector<VectorXd> sigma_points = generateSigmaPoints(aug_state, aug_cov);
            std::vector<double> sigma_weights = generateSigmaWeights(num_aug);

            std::vector<VectorXd> z_sig;
            for(const auto& sigma_point : sigma_points){z_sig.push_back(lidarMeasurementModel(sigma_point, map_beacon.x, map_beacon.y));}

            VectorXd z_mean = VectorXd::Zero(num_w);
            for(int i = 0; i<z_sig.size();++i)
            {
                z_mean += sigma_weights[i]*z_sig[i];
            }

            MatrixXd S = MatrixXd::Zero(num_w, num_w);
            for(unsigned int i = 0; i<z_sig.size(); ++i)
            {
                VectorXd z_diff = normaliseLidarMeasurement(z_sig[i] - z_mean);

                S += sigma_weights[i] * z_diff * z_diff.transpose();
            }

            MatrixXd cross_cov = MatrixXd::Zero(num_state, num_w);

            for(unsigned int i = 0; i<2*num_state+1;++i)
            {
                VectorXd x_diff = normaliseState(sigma_points[i].head(num_state) - state);
                VectorXd z_diff = normaliseLidarMeasurement(z_sig[i] - z_mean);
                
                cross_cov += sigma_weights[i] * x_diff * z_diff.transpose();
            }

            VectorXd y = normaliseLidarMeasurement(z - z_mean);
            MatrixXd K = cross_cov * S.inverse();

            state = state + K * y;
            cov = cov - K * S * K.transpose();
            
           

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
        // Hint: You can use the constants: ACCEL_STD, GYRO_STD
        // HINT: use the wrapAngle() function on angular values to always keep angle
        // values within correct range, otherwise strange angle effects might be seen.
        // ----------------------------------------------------------------------- //
        // ENTER YOUR CODE HERE
        size_t n_x = state.size();
        int n_w = 2;
        int n_aug = n_x + n_w;
        VectorXd aug_state = VectorXd::Zero(n_aug);
        aug_state.head(n_x) = state;

        MatrixXd Q = MatrixXd::Zero(n_w, n_w);
        Q(0,0) = GYRO_STD*GYRO_STD;
        Q(1,1) = ACCEL_STD*ACCEL_STD;

        MatrixXd aug_cov = MatrixXd::Zero(n_aug, n_aug);
        aug_cov.topLeftCorner(n_x, n_x) = cov;
        aug_cov.bottomRightCorner(n_w, n_w) = Q;

        std::vector<VectorXd> sigmaPoints = generateSigmaPoints(aug_state, aug_cov);
        std::vector<double> sigmaWeights = generateSigmaWeights(n_aug);

        std::vector<VectorXd> predicted_sigma_points;

        for(const auto& sigmaPoint: sigmaPoints)
        {
            
            predicted_sigma_points.push_back(vehicleProcessModel(sigmaPoint, gyro.psi_dot, dt));
        }

        state = VectorXd::Zero(n_x);
        for(int i=0; i<predicted_sigma_points.size();++i)
        {
            state += sigmaWeights[i]*predicted_sigma_points[i];
        }
        state = normaliseState(state);

        cov = MatrixXd::Zero(n_x, n_x);
        for(int i=0; i<predicted_sigma_points.size();++i)
        {
            VectorXd D = normaliseState(predicted_sigma_points[i]-state);
            cov += sigmaWeights[i]* D * D.transpose();
        }


        // ----------------------------------------------------------------------- //

        setState(state);
        setCovariance(cov);
    } 
}

void KalmanFilter::handleGPSMeasurement(GPSMeasurement meas)
{
    if(isInitialised())
    {
        VectorXd state = getState();
        MatrixXd cov = getCovariance();

        VectorXd z = Vector2d::Zero();
        MatrixXd H = MatrixXd(2,4);
        MatrixXd R = Matrix2d::Zero();

        z << meas.x,meas.y;
        H << 1,0,0,0,0,1,0,0;
        R(0,0) = GPS_POS_STD*GPS_POS_STD;
        R(1,1) = GPS_POS_STD*GPS_POS_STD;

        VectorXd z_hat = H * state;
        VectorXd y = z - z_hat;
        MatrixXd S = H * cov * H.transpose() + R;
        MatrixXd K = cov*H.transpose()*S.inverse();

        state = state + K*y;
        cov = (Matrix4d::Identity() - K*H) * cov;

        setState(state);
        setCovariance(cov);
    }
    else
    {
        // You may modify this initialisation routine if you can think of a more
        // robust and accuracy way of initialising the filter.
        // ----------------------------------------------------------------------- //
        // YOU ARE FREE TO MODIFY THE FOLLOWING CODE HERE

        VectorXd state = Vector4d::Zero();
        MatrixXd cov = Matrix4d::Zero();

        state(0) = meas.x;
        state(1) = meas.y;
        cov(0,0) = GPS_POS_STD*GPS_POS_STD;
        cov(1,1) = GPS_POS_STD*GPS_POS_STD;
        cov(2,2) = INIT_PSI_STD*INIT_PSI_STD;
        cov(3,3) = INIT_VEL_STD*INIT_VEL_STD;

        setState(state);
        setCovariance(cov);

        // ----------------------------------------------------------------------- //
    }             
}

void KalmanFilter::handleLidarMeasurements(const std::vector<LidarMeasurement>& dataset, const BeaconMap& map)
{
    // Assume No Correlation between the Measurements and Update Sequentially
    for(const auto& meas : dataset) {handleLidarMeasurement(meas, map);}
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
