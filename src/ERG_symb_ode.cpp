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

//vector_type x(14, 0.0); //Initial condition
double q_des_init[7] = {0,0,0,0,0,0,0};
double T_sim = 100;
int n_iter = 0;
MatrixXd q_out_vec = MatrixXd::Zero(5000,28);
int my_flag = 0;



class KUKA_INVKIN {
	public:
		KUKA_INVKIN();
		void run();
		bool init_robot_model();
		void ctrl_loop();
                VectorXd estimated_dynamic(VectorXd x,VectorXd q_des);
                VectorXd dynamic_mdl(VectorXd x);
                MatrixXd num_J(VectorXd x);
                VectorXd NF_att(VectorXd r_,VectorXd v_);
                VectorXd NF_rep(VectorXd r_);
                double DSM(VectorXd q_in,VectorXd _v);
                double eval_Delta(VectorXd x,VectorXd lim);
                MatrixXd Kp = MatrixXd::Zero(7, 7);
	        MatrixXd Kd = MatrixXd::Zero(7, 7);
                void joint_states_cb( sensor_msgs::JointState );
                
                

	private:
      	   	
		ros::NodeHandle _nh;
                VectorXd qlim;
                VectorXd qdlim;
                VectorXd ulim;
 
                bool _first_js;
                VectorXd _q_in = VectorXd::Zero(7);
                VectorXd _dq_in = VectorXd::Zero(7);
                ros::Subscriber _my_sub;
                ros::Publisher _my_pub_v;

               

};



bool KUKA_INVKIN::init_robot_model() {
	//std::string robot_desc_string;
	//_nh.param("robot_description", robot_desc_string, std::string());


        Kp(0, 0) = 30+100*0;
        Kp(1, 1) = 30+100*0;
        Kp(2, 2) = 30+100*0;
        Kp(3, 3) = 30+100*0;
        Kp(4, 4) = 10+100*0;
        Kp(5, 5) = 10+100*0;
        Kp(6, 6) = 10+100*0;

        Kd(0, 0) = 1+10*0;
        Kd(1, 1) = 1+10*0;
        Kd(2, 2) = 1+10*0;
        Kd(3, 3) = 1+10*0;
        Kd(4, 4) = 0.1+10*0;
        Kd(5, 5) = 0.1+10*0;
        Kd(6, 6) = 0.1+10*0;

        qlim = VectorXd::Zero(7);
        qdlim = VectorXd::Zero(7);
        ulim = VectorXd::Zero(7);
        
        qlim << 170*3.14/180,120*3.14/180,170*3.14/180,120*3.14/180,170*3.14/180,120*3.14/180,175*3.14/180;
        qdlim << 85*3.14/180,85*3.14/180,100*3.14/180,75*3.14/180,130*3.14/180,135*3.14/180,135*3.14/180;
        ulim << 320,320,170,170,110,40,40;


	return true;
}



KUKA_INVKIN::KUKA_INVKIN() {

	if (!init_robot_model()) exit(1);	
        
        _my_pub_v = _nh.advertise< std_msgs::Float64MultiArray > ("/ref",0);  
        _my_sub = _nh.subscribe("/lbr_iiwa/joint_states", 0, &KUKA_INVKIN::joint_states_cb, this);

        _first_js = false;

}


void KUKA_INVKIN::joint_states_cb( sensor_msgs::JointState js ) {

	for(int i=0; i<7; i++ ) { 
		_q_in(i)= js.position[i];
		_dq_in(i) = js.velocity[i];
		
	}
        //cout<<_q_in.transpose()<<endl;
	_first_js = true;

}

VectorXd KUKA_INVKIN::NF_att(VectorXd r_,VectorXd v_){

/*ERG Parameters*/

VectorXd ro_att = VectorXd::Zero(7); 
double ro1_att = 0;
VectorXd ro2_att = VectorXd::Zero(2);
double eta = 0.00001; 
VectorXd ro3_att = VectorXd::Zero(7);


/**/

ro_att = r_ - v_;
ro1_att = ro_att.norm();
cout<<"norm r-v: "<<ro1_att<<endl;
ro2_att(0) = ro1_att;
ro2_att(1) = eta; 
ro3_att = ro_att/ro2_att.maxCoeff(); 

return ro3_att;

}

VectorXd KUKA_INVKIN::NF_rep(VectorXd r_){

VectorXd ro_rep = VectorXd::Zero(7);

VectorXd Temp1 = VectorXd::Zero(2);
VectorXd Temp2 = VectorXd::Zero(2);

double zeta = 0.15;
double delta = 0.1;
double temp1 = 0;
double temp2 = 0;


for(int i=0;i<7;i++){

temp1 = abs(r_(i)+qlim(i));
temp1 = (zeta - temp1)/(zeta - delta);
Temp1(0) = temp1;
Temp1(1) = 0;

temp2 = abs(r_(i)-qlim(i));
temp2 = (zeta - temp2)/(zeta - delta);
Temp2(0) = temp2;
Temp2(1) = 0;

ro_rep(i) = Temp1.maxCoeff() - Temp2.maxCoeff();

}

return ro_rep;

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


		for (int i = 0; i < 14; i++) {
                          x2(i) = x[i];
                          q_out_vec(my_flag,i) = x2(i);
                }
		J_t = ik.num_J(x2);
		//cout << "J: " << J_t << endl;
		for (int i = 0; i < 14; i++) {
			for (int j = 0; j < 14; j++) {
				J(i, j) = J_t(i, j);
			}
		}
		
		for (int i = 0; i < 14; i++) dfdt[i] = 0.0;
                my_flag = my_flag + 1;
		
};


VectorXd KUKA_INVKIN::estimated_dynamic(VectorXd q_in,VectorXd q_des_){

        vector_type x(14, 0.0); //Initial condition
        VectorXd q_out = VectorXd::Zero(14);

        for(int i =0;i<14;i++) x[i] = q_in(i);
        //double my_time = 0;

        for (int i = 0; i < 7; i++) q_des_init[i] = q_des_(i);

        //cout<<"estimate_dynamic"<<endl;

        //auto start = high_resolution_clock::now();
        //cout<<"start_time: "<<my_time<<endl;
	size_t num_of_steps = integrate_adaptive(make_dense_output< rosenbrock4< double > >(1.0e-3, 1.0e-6),
		make_pair(stiff_system,stiff_system_jacobi),
		x, 0.0, T_sim, 0.001);//,
                //cout << phoenix::arg_names::arg2 << " " << phoenix::arg_names::arg1[0] << "\n" );
               
        //auto stop = high_resolution_clock::now();                
        //auto duration = duration_cast<microseconds>(stop - start);
        //my_time = duration.count()*0.000001;
        //cout<<"time: "<<my_time<<endl;

        int my_flag = 0;
        n_iter = num_of_steps;


        for(int i =0;i<14;i++) q_out(i) = x[i];
        return q_out;

}


double KUKA_INVKIN::DSM(VectorXd q_in,VectorXd _v){

/*INITIALIZATION*/
VectorXd _q = VectorXd::Zero(7);
VectorXd _dq = VectorXd::Zero(7);

for(int i =0;i<7;i++) _q(i) = q_in(i);
for(int i =0;i<7;i++) _dq(i) = q_in(i+7);




// ERG Tuning Parameters /*TO TUNE*/
double kq = 8;
double kqd = 5;
double ktao = 0.03;
double k_e = 85;
double E_term = 150;


//Vector
VectorXd q_out = VectorXd::Zero(14);
VectorXd tao = VectorXd::Zero(7);
VectorXd e = VectorXd::Zero(7);
VectorXd out = VectorXd::Zero(14); //out = [q_ dq_]
VectorXd q_ = VectorXd::Zero(7);
VectorXd dq_ = VectorXd::Zero(7);


VectorXd Delta_vec = VectorXd::Zero(3);
VectorXd Delta_vec2 = VectorXd::Zero(2);
VectorXd Delta_vec3 = VectorXd::Zero(2);


double Delta1 = 0;
VectorXd V = VectorXd::Zero(1);
double Delta_V = 0;
double Delta = 0;



//For dynamic prediction
VectorXd G = VectorXd::Zero(7);
MatrixXd B = MatrixXd::Zero(7,7);

//before we have to set the reference

q_out = estimated_dynamic(q_in,_v);
//cout<<"n_iter: "<<n_iter<<endl;
//Vector
VectorXd Delta_q = VectorXd::Zero(n_iter);
VectorXd Delta_dq = VectorXd::Zero(n_iter);
VectorXd Delta_tao = VectorXd::Zero(n_iter);

for(int i=0;i<n_iter;i++){

q_out = q_out_vec.row(i);
for (int j=0;j<7;j++) q_(j) = q_out(j); //unpack position
for (int j=0;j<7;j++) dq_(j) = q_out(j+7); //unpack velocity

G = my_matrix_G(_v);
//torque 
e = _v - q_; //error in position 
tao = -Kd*dq_ + Kp*e + G;

Delta_q(i) = eval_Delta(q_,qlim);
Delta_dq(i) = eval_Delta(dq_,qdlim);
Delta_tao(i) = eval_Delta(tao,ulim);

}

Delta_vec(0) =  kq*Delta_q.minCoeff();
Delta_vec(1) =  kqd*Delta_dq.minCoeff(); 
Delta_vec(2) =  ktao*Delta_tao.minCoeff();   


Delta1 = Delta_vec.minCoeff();

/*Lyapunov Function evautation*/
B = my_matrix_B(q_);

V = 0.5*dq_.transpose()*B*dq_ + 0.5*e.transpose()*Kp*e;
Delta_V = k_e*(E_term - V(0));  


Delta_vec2(0) = Delta1;
Delta_vec2(1) = Delta_V;

Delta_vec3(0) = Delta_vec2.minCoeff();
Delta_vec3(1) = 0;

Delta = Delta_vec3.maxCoeff();


return Delta;

}


double KUKA_INVKIN::eval_Delta(VectorXd x,VectorXd lim){

VectorXd Temp = VectorXd::Zero(2);
VectorXd Delta_v = VectorXd::Zero(7);
double Delta_out = 0;

for(int k = 0; k<7;k++){
Temp(0) = x(k) + lim(k);
Temp(1) = - x(k) + lim(k);
Delta_v(k) = Temp.minCoeff();
}

Delta_out = Delta_v.minCoeff();
}





void KUKA_INVKIN::ctrl_loop() {

/*vector inti*/
VectorXd q_in = VectorXd::Zero(14);
VectorXd q_out = VectorXd::Zero(14);
VectorXd r_ = VectorXd::Zero(7);
VectorXd v_ = VectorXd::Zero(7);



/*ERG Parameters*/
VectorXd ro_att = VectorXd::Zero(7); 
VectorXd ro_rep = VectorXd::Zero(7); 
VectorXd ro = VectorXd::Zero(7);
std_msgs::Float64MultiArray cmd_pub_v;
cmd_pub_v.data.resize(7);

float Ts = 0.001;

double deg2rad = 3.14/180;
float offset = 50;
double Delta = 0; 

while(!_first_js) usleep(0.1);
_first_js = false;
for(int i=0; i<7; i++ ) q_in(i) = _q_in(i);
for(int i=0; i<7; i++ ) q_in(i+7) = _dq_in(i);


for(int i=0; i<7; i++ ) r_(i)  = _q_in(i)+offset*deg2rad;

for(int i = 0;i<7;i++) v_(i) = q_in(i); //initialize the v at the same q0


int flag = 0;

double my_time = 0;


cout<<"initial position: "<<_q_in.transpose()<<endl;
cout<<"initial input: "<<v_.transpose()<<endl;
cout<<"reference: "<<r_.transpose()<<endl;

cin>>flag;

ros::Rate r(100);
while( ros::ok() ) {

while(!_first_js)
_first_js = false;
for(int i=0; i<7; i++ ) q_in(i) = _q_in(i);
for(int i=0; i<7; i++ ) q_in(i+7) = _dq_in(i);



ro_att = KUKA_INVKIN::NF_att(r_,v_); //attraction field
//cout<<"ro_att: "<<ro_att.transpose()<<endl; 
ro_rep = KUKA_INVKIN::NF_rep(r_); //repulsion field
//cout<<"ro_rep: "<<ro_rep.transpose()<<endl;
ro = ro_att + ro_rep; //final navigation field

T_sim = 100*Ts;
auto start = high_resolution_clock::now();
Delta = KUKA_INVKIN::DSM(q_in,v_); //Dynamic Safety Margin Function 
auto stop = high_resolution_clock::now();                
auto duration = duration_cast<microseconds>(stop - start);
my_time = duration.count()*0.000001;
cout<<"time: "<<my_time<<endl;

v_ = ro*Delta*Ts + v_;

cout<<"applied v: "<<v_.transpose()<<endl;
for(int i=0; i<7; i++) cmd_pub_v.data[i] = v_(i);
_my_pub_v.publish(cmd_pub_v);


r.sleep(); 
}



}


//T_sim = Ts;
//q_out= estimated_dynamic(q_in,v_);
//q_in = q_out;


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
