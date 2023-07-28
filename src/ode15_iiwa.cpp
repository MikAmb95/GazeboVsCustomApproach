#include "ros/ros.h"
#include "boost/thread.hpp"
#include "sensor_msgs/JointState.h"
#include "geometry_msgs/Pose.h"
#include <std_msgs/Float64.h>
#include <std_msgs/Float64MultiArray.h>
#include <iostream>
#include <sstream>
#include <fstream>
#include <stdio.h>
#include <fcntl.h>
#include <arpa/inet.h>
#include <netinet/in.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <sys/mman.h>
#include <sys/io.h>
#include <sys/time.h>
#include <netdb.h>

#include <chrono>



#include <ros/package.h>
#include <tf/tf.h>
#include <tf_conversions/tf_eigen.h>

#include <Eigen/Dense>

#include "my_roscpp_library/my_roscpp_library.h"

#include <boost/numeric/odeint.hpp>
#include <boost/phoenix/core.hpp>

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>


using namespace Eigen;
using namespace std;
using namespace std::chrono;

using namespace std::chrono;
using namespace boost::numeric::odeint;
namespace phoenix = boost::phoenix;

typedef std::vector< double > state_type;
//[ stiff_system_definition
typedef boost::numeric::ublas::vector< double > vector_type;
typedef boost::numeric::ublas::matrix< double > matrix_type;

vector_type x(14, 0.0); //Initial condition
double q_des_init[7] = {0,0,0,0,0,0,0};






class KUKA_INVKIN {
	public:
		KUKA_INVKIN();
		void run();
		bool init_robot_model();
		void ctrl_loop();
                void estimated_dynamic();
                VectorXd dynamic_mdl(VectorXd x);
                MatrixXd num_J(VectorXd x);
                MatrixXd Kp = MatrixXd::Zero(7, 7);
	        MatrixXd Kd = MatrixXd::Zero(7, 7);
                
                

	private:
      	   	
		ros::NodeHandle _nh;
               

};



bool KUKA_INVKIN::init_robot_model() {
	//std::string robot_desc_string;
	//_nh.param("robot_description", robot_desc_string, std::string());

	Kp(0, 0) = 100;
	Kp(1, 1) = 100;
	Kp(2, 2) = 100;
	Kp(3, 3) = 100;
	Kp(4, 4) = 100;
	Kp(5, 5) = 100;
	Kp(6, 6) = 100;

	Kd(0, 0) = 10;
	Kd(1, 1) = 10;
	Kd(2, 2) = 10;
	Kd(3, 3) = 10;
	Kd(4, 4) = 10;
	Kd(5, 5) = 10;
	Kd(6, 6) = 10;
	return true;
}



KUKA_INVKIN::KUKA_INVKIN() {

	if (!init_robot_model()) exit(1);	
}




VectorXd KUKA_INVKIN::dynamic_mdl(VectorXd x) {
        


	VectorXd q = VectorXd::Zero(7);
	VectorXd dq = VectorXd::Zero(7);
	VectorXd ex = VectorXd::Zero(14);
        MatrixXd B = MatrixXd::Zero(7, 7);
	MatrixXd C = MatrixXd::Zero(7, 7);
	MatrixXd G = VectorXd::Zero(7);
        VectorXd tau = VectorXd::Zero(7);
	VectorXd ddq = VectorXd::Zero(7);
        VectorXd q_des = VectorXd::Zero(7);


	for (int i = 0; i < 7; i++)q(i) = x(i);
	for (int i = 0; i < 7; i++)dq(i) = x(i + 7);
        for (int i = 0; i < 7; i++) q_des(i) = q_des_init[i];



	B = my_matrix_B(q);
        C = my_matrix_C(q, dq);
        G = my_matrix_G(q);
        
        tau = Kp * (q_des - q) - Kd * dq + G;
        ddq = B.inverse() * (-C * dq - G + tau);

        
	for (int i = 0; i < 7; i++) ex(i) = x(i+7);

	for (int i = 0; i < 7; i++) ex(i+7) = ddq(i);
        
	return ex;


}

MatrixXd KUKA_INVKIN::num_J(VectorXd x)
{

        int dim_s = 14;
	double eps = 0.0000001;
	MatrixXd J2 = MatrixXd::Zero(14, 14);
	VectorXd fx_p = VectorXd::Zero(14);
	VectorXd fx = dynamic_mdl(x);
        KUKA_INVKIN ik;
	VectorXd x_p = x;

	for (int i = 0; i < dim_s; i++) {
		x_p(i) = x_p(i) + eps;
		fx_p = dynamic_mdl(x_p);
		for (int j = 0; j < 14; j++) {
			J2(j, i) = (fx_p(j) - fx(j)) / eps;
		}
		x_p(i) = x(i);
	}
	return J2;

}


void stiff_system(const vector_type& x, vector_type& dxdt, double /* t */)
	{


                KUKA_INVKIN ik;
		MatrixXd B = MatrixXd::Zero(7, 7);
		MatrixXd C = MatrixXd::Zero(7, 7);
		MatrixXd G = VectorXd::Zero(7);
		VectorXd q = VectorXd::Zero(7);
		VectorXd dq = VectorXd::Zero(7);
                VectorXd ddq = VectorXd::Zero(7);
		VectorXd tau = VectorXd::Zero(7);
                VectorXd q_des = VectorXd::Zero(7);

                for (int i = 0; i < 7; i++) q_des(i) = q_des_init[i];
                for (int i = 0; i < 7; i++) q(i) = x[i];
                for (int i = 0; i < 7; i++) dq(i) = x[i+7];

		B = my_matrix_B(q);
		C = my_matrix_C(q, dq);
		G = my_matrix_G(q);
		tau = ik.Kp * (q_des - q) - ik.Kd * dq + G;
		ddq = B.inverse() * (-C * dq - G + tau);
                
                for (int i = 0; i < 7; i++) dxdt[i] = x[i+7];
                for (int i = 0; i < 7; i++) dxdt[i+7] = ddq(i);

};

void stiff_system_jacobi(const vector_type& x, matrix_type& J, const double& /* t */, vector_type& dfdt)
	{
		MatrixXd J_t = MatrixXd::Zero(14, 14);
		VectorXd x2 = VectorXd::Zero(14);
                KUKA_INVKIN ik;


		for (int i = 0; i < 14; i++)x2(i) = x[i];

		J_t = ik.num_J(x2);
		//cout << "J: " << J_t << endl;
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 14; j++) {
				J(i, j) = J_t(i, j);
			}
		}
		
		for (int i = 0; i < 14; i++) dfdt[i] = 0.0;
		
};


void KUKA_INVKIN::estimated_dynamic(){


        double my_time = 0;

        auto start = high_resolution_clock::now();
        //cout<<"start_time: "<<my_time<<endl;
	size_t num_of_steps = integrate_adaptive(make_dense_output< rosenbrock4< double > >(1.0e-3, 1.0e-6),
		make_pair(stiff_system,stiff_system_jacobi),
		x, 0.0, 100.0, 0.001);//,
                //cout << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << "\n" );
               
        auto stop = high_resolution_clock::now();                
        auto duration = duration_cast<microseconds>(stop - start);
        my_time = duration.count()*0.000001;
        cout<<"time: "<<my_time<<endl;

}


void KUKA_INVKIN::ctrl_loop() {

        //ros::Rate r(1);

       
        

        //while( ros::ok() ) {
        VectorXd q_des = VectorXd::Zero(7);
        VectorXd q_0 = VectorXd::Zero(14);
        VectorXd dq_0 = VectorXd::Zero(7);



        int flag = 0;
        
        cout<<"first dynamic_simulation"<<endl;
        for(int i=0;i<14;i++) x[i]=q_0(i);
        cout<<"Initial q: ";
        for(int i=0;i<14;i++) cout<<x[i]<<" ";
        cout<<endl;
        q_des<<10*3.14/180,5*3.14/180,-5*3.14/180,20*3.14/180,-10*3.14/180,5*3.14/180,-5*3.14/180;
        for (int i = 0; i < 7; i++) q_des_init[i] = q_des(i);
        estimated_dynamic();
        cout<<"end simulation"<<endl;
        cout<<"Final q: ";
        for(int i=0;i<14;i++) cout<<x[i]<<" ";
        cout<<endl;
        cout<<"qd: ";
        for(int i=0;i<7;i++) cout<<q_des_init[i]<<" ";
        cout<<endl;
        cin>>flag;
        
        cout<<"second dynamic_simulation"<<endl;
        //for(int i=0;i<14;i++) q_0(i) = x[i];
        //for(int i=0;i<14;i++) x[i] = q_0(i);
        cout<<"Initial q: ";
        for(int i=0;i<14;i++) cout<<x[i]<<" ";
        cout<<endl;
        q_des<<30*3.14/180,-5*3.14/180,-15*3.14/180,25*3.14/180,-5*3.14/180,-5*3.14/180,10*3.14/180;
        for (int i = 0; i < 7; i++) q_des_init[i] = q_des(i);
        estimated_dynamic();
        cout<<"end simulation"<<endl;
        cout<<"Final q: ";
        for(int i=0;i<14;i++) cout<<x[i]<<" ";
        cout<<endl;
        cout<<"qd: ";
        for(int i=0;i<7;i++) cout<<q_des_init[i]<<" ";
        cout<<endl;
        cin>>flag;
       


        //r.sleep();
        //}

}


void KUKA_INVKIN::run() {

	boost::thread ctrl_loop_t ( &KUKA_INVKIN::ctrl_loop, this);
	ros::spin();	

}




int main(int argc, char** argv) {
         
        cout<<"StartOde15"<<endl;
	ros::init(argc, argv, "ode15");
	KUKA_INVKIN ik;
	ik.run();

	return 0;
}
