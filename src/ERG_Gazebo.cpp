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

#include <kdl_parser/kdl_parser.hpp>
#include <kdl/chainfksolvervel_recursive.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolverpos_nr.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainjnttojacdotsolver.hpp>
#include <kdl/chaindynparam.hpp>
#include <ros/package.h>
#include <tf/tf.h>
#include <tf_conversions/tf_eigen.h>
#include <Eigen/Dense>


using namespace Eigen;
using namespace std;
using namespace KDL;

using namespace std::chrono;

class KUKA_INVKIN {
	public:
		KUKA_INVKIN();
		void run();
		bool init_robot_model();
                void pd_g();
		void ctrl_loop();
                double ERG_function(VectorXd _q,VectorXd _dq,VectorXd _v);
                VectorXd forw_dy2(VectorXd _q,VectorXd _dq,VectorXd _v);
                VectorXd NF_att(VectorXd r_,VectorXd v_);
                VectorXd NF_rep(VectorXd r_);
                double eval_Delta(VectorXd x,VectorXd lim);
                void joint_states_cb( sensor_msgs::JointState );


	private:
	
		ros::NodeHandle _nh;
		KDL::Tree iiwa_tree;
		KDL::Chain _k_chain;

                VectorXd _q_in = VectorXd::Zero(7);
                VectorXd _dq_in = VectorXd::Zero(7);

                bool _first_js;

                MatrixXd Kd; 	
		MatrixXd Kp;
                KDL::ChainDynParam *_dyn_param;
                ros::Publisher _my_pub_q;
                ros::Publisher _my_pub_v;
                ros::Subscriber _my_sub;

                VectorXd qlim;
                VectorXd qdlim;
                VectorXd ulim;

               boost::shared_ptr<sensor_msgs::JointState const>msg_r;

                

	
};


bool KUKA_INVKIN::init_robot_model() {
	std::string robot_desc_string;
	_nh.param("robot_description", robot_desc_string, std::string());
	if (!kdl_parser::treeFromString(robot_desc_string, iiwa_tree)){
		ROS_ERROR("Failed to construct kdl tree");
		return false;
	}

	std::string base_link = "lbr_iiwa_link_0";
	std::string tip_link  = "lbr_iiwa_link_7";
	if ( !iiwa_tree.getChain(base_link, tip_link, _k_chain) ) return false;

        _dyn_param = new KDL::ChainDynParam(_k_chain,KDL::Vector(0,0,-9.81));

        Kp = MatrixXd::Zero(7,7);
        Kd = MatrixXd::Zero(7,7);
        Kp(0, 0) = 30;
        Kp(1, 1) = 30;
        Kp(2, 2) = 30;
        Kp(3, 3) = 30;
        Kp(4, 4) = 10;
        Kp(5, 5) = 10;
        Kp(6, 6) = 10;

        Kd(0, 0) = 1;
        Kd(1, 1) = 1;
        Kd(2, 2) = 1;
        Kd(3, 3) = 1;
        Kd(4, 4) = 0.1;
        Kd(5, 5) = 0.1;
        Kd(6, 6) = 0.1;

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
	ROS_INFO("Robot tree correctly loaded from parameter server!");

	cout << "Joints and Link: " << iiwa_tree.getNrOfJoints() << " - " << iiwa_tree.getNrOfSegments() << endl;
        _my_pub_q = _nh.advertise< std_msgs::Float64MultiArray > ("/iiwa_q",0); 
        _my_pub_v = _nh.advertise< std_msgs::Float64MultiArray > ("/ref",0);  
        _my_sub = _nh.subscribe("/lbr_iiwa/joint_states", 0, &KUKA_INVKIN::joint_states_cb, this);

        _first_js = false;


}


void KUKA_INVKIN::joint_states_cb( sensor_msgs::JointState js ) {

	for(int i=0; i<7; i++ ) { 
		_q_in(i)= js.position[i];
		_dq_in(i) = js.velocity[i];
		
	}
	_first_js = true;

}

VectorXd KUKA_INVKIN::forw_dy2(VectorXd _q,VectorXd _dq,VectorXd _v){

    KDL::JntArray coriol_(7);
    KDL::JntArray grav_(7);
    KDL::JntSpaceInertiaMatrix jsim_;
    jsim_.resize(_k_chain.getNrOfJoints());

    KDL::JntArray *q;
    KDL::JntArray *dq;

    q = new KDL::JntArray( _k_chain.getNrOfJoints() );
    dq = new KDL::JntArray( _k_chain.getNrOfJoints() );

    VectorXd out = VectorXd::Zero(14);
    VectorXd _ddq = VectorXd::Zero(7);
    VectorXd q_tilda = VectorXd::Zero(7);
    VectorXd tao = VectorXd::Zero(7);
    double Ts = 0.001;

    MatrixXd m = MatrixXd::Zero(7,7);



    for(int i = 0;i<7;i++){
    q->data[i] = _q(i);
    dq->data[i] = _dq(i);
    }

    _dyn_param->JntToMass(*q, jsim_);
    _dyn_param->JntToCoriolis(*q, *dq, coriol_);
    _dyn_param->JntToGravity(*q, grav_);

    m = jsim_.data;
    m = m.inverse();    

    q_tilda = _v - _q;
    tao =  Kp*q_tilda - Kd*_dq + grav_.data;
    _ddq = m*(- coriol_.data - grav_.data + tao);

    _dq = Ts*_ddq + _dq;
    _q = Ts*_dq + _q;

    for (int i=0;i<7;i++) out(i) = _q(i);
    for (int i=0;i<7;i++) out(i+7) = _dq(i);

    return out;
                       
}

double KUKA_INVKIN::ERG_function(VectorXd _q,VectorXd _dq,VectorXd _v){


int flag = 0;

// ERG Tuning Parameters /*TO TUNE*/
double kq = 8;
double kqd = 5;
double ktao = 0.03;
double k_e = 85;
double E_term = 150;
int iter = 100*0+10000; //number of interation for prediction 

//Vector
VectorXd Delta_q = VectorXd::Zero(iter);
VectorXd Delta_dq = VectorXd::Zero(iter);
VectorXd Delta_tao = VectorXd::Zero(iter);

VectorXd Delta_vec = VectorXd::Zero(3);
VectorXd Delta_vec2 = VectorXd::Zero(2);
VectorXd Delta_vec3 = VectorXd::Zero(2);


double Delta1 = 0;
VectorXd V = VectorXd::Zero(1);
double Delta_V = 0;
double Delta = 0;


VectorXd tao = VectorXd::Zero(iter);
VectorXd e = VectorXd::Zero(7);


 	
VectorXd out = VectorXd::Zero(14); //out = [q_ dq_]
VectorXd q_ = VectorXd::Zero(7);
VectorXd dq_ = VectorXd::Zero(7);




//For dynamic prediction
KDL::JntArray grav_(7);
KDL::JntArray *q_in;
q_in = new KDL::JntArray( _k_chain.getNrOfJoints() );

KDL::JntSpaceInertiaMatrix jsim_;
jsim_.resize(_k_chain.getNrOfJoints());
MatrixXd B = MatrixXd::Zero(7,7);;

/*-------------*/


for(int i = 0;i<iter;i++){
out = KUKA_INVKIN::forw_dy2(_q,_dq,_v); //predict dynamic for a inter-steps out = [q_ dq_]


for (int j=0;j<7;j++) q_(j) = out(j); //unpack position
for (int j=0;j<7;j++) dq_(j) = out(j+7); //unpack velocity


//gravity compensation
for(int k = 0;k<7;k++) q_in->data[k] = _v(k); //to evalueta the gravity comp
_dyn_param->JntToGravity(*q_in, grav_);

//torque 
e = _v - _q; //error in position 
tao = -Kd*_dq + Kp*e + grav_.data;

Delta_q(i) = eval_Delta(_q,qlim);
Delta_dq(i) = eval_Delta(_dq,qdlim);
Delta_tao(i) = eval_Delta(tao,ulim);

_q = q_; //re-assing position for next iteration
_dq = dq_; //re-assign velocity for next iteration
}


Delta_vec(0) =  kq*Delta_q.minCoeff();
Delta_vec(1) =  kqd*Delta_dq.minCoeff(); 
Delta_vec(2) =  ktao*Delta_tao.minCoeff();   


Delta1 = Delta_vec.minCoeff();

/*Lyapunov Function evautation*/
for(int k = 0;k<7;k++) q_in->data[k] = _q(k);
_dyn_param->JntToMass(*q_in, jsim_);
B = jsim_.data;

V = 0.5*_dq.transpose()*B*_dq + 0.5*e.transpose()*Kp*e;
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
cout<<"norm r-v"<<ro1_att<<endl;
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

void KUKA_INVKIN::ctrl_loop() {

VectorXd out = VectorXd::Zero(14);
VectorXd q_ = VectorXd::Zero(7);
VectorXd dq_ = VectorXd::Zero(7);
VectorXd v_ = VectorXd::Zero(7);
VectorXd r_ = VectorXd::Zero(7);

/*ERG Parameters*/
VectorXd ro_att = VectorXd::Zero(7); 
VectorXd ro_rep = VectorXd::Zero(7); 
VectorXd ro = VectorXd::Zero(7);
double Delta = 0; 
/**/

std_msgs::Float64MultiArray cmd_pub_q;
cmd_pub_q.data.resize(7);

std_msgs::Float64MultiArray cmd_pub_v;
cmd_pub_v.data.resize(7);

float Ts = 0.001;

double deg2rad = 3.14/180;
float offset = 50;



//r_ <<0,0.174,0,-1.6,0,1.4,0.174;
q_ <<0,0.172,0,-1.4,0,1.5,0;
//msg_r = ros::topic::waitForMessage<sensor_msgs::JointState>("/lbr_iiwa/joint_states"); //to know the robot status
//for(int i=0; i<7; i++ ) q_(i) = msg_r->position[i];

int flag;
double time_ERG = 0;

//while(!_first_js) usleep(0.1);
//_first_js = false;
//q_ = _q_in;
//dq_ = _dq_in;



for(int i=0; i<7; i++ ) r_(i)  = q_(i)+offset*deg2rad;

cout<<"desiderd reference: "<<r_.transpose()<<endl;
v_ = q_; //initilialize the applied reference at the same initial condition
//r_ = v_;
       ros::Rate r(50);
       while( ros::ok() ) {	
        auto start = high_resolution_clock::now(); 
        //while(!_first_js)
        //_first_js = false;
        //q_ = _q_in;
        //dq_ = _dq_in;

        ro_att = KUKA_INVKIN::NF_att(r_,v_);
        //cout<<"ro_att: "<<ro_att.transpose()<<endl; 
        ro_rep = KUKA_INVKIN::NF_rep(r_);
        //cout<<"ro_rep: "<<ro_rep.transpose()<<endl;
        ro = ro_att + ro_rep;


        //auto start = high_resolution_clock::now(); 


        Delta = KUKA_INVKIN::ERG_function(q_,dq_,v_);


        //cout<<"Delta: "<<Delta<<endl;
        v_ = ro*Delta*Ts + v_;
        
        //cout<<"applied v: "<<v_.transpose()<<endl;
        //for(int i=0; i<7; i++) cmd_pub_v.data[i] = v_(i);
        //_my_pub_v.publish(cmd_pub_v);

        auto stop = high_resolution_clock::now();
        auto duration = duration_cast<microseconds>(stop - start);
        time_ERG = duration.count()*0.000001;
        cout<<"time: "<<time_ERG<<endl;

        out = KUKA_INVKIN::forw_dy2(q_,dq_,v_); //predict dynamic for a inter-steps out = [q_ dq_]
        for (int j=0;j<7;j++) q_(j) = out(j); //unpack position
        for (int j=0;j<7;j++) dq_(j) = out(j+7); //unpack velocity

        r.sleep();             
	}
}





void KUKA_INVKIN::run() {

	boost::thread ctrl_loop_t ( &KUKA_INVKIN::ctrl_loop, this);
	ros::spin();	

}




int main(int argc, char** argv) {

        cout<<"ERG"<<endl;
	ros::init(argc, argv, "iiwa_kdl_erg");
	KUKA_INVKIN ik;
	ik.run();

	return 0;
}


 /*
//msg_r = ros::topic::waitForMessage<sensor_msgs::JointState>("/lbr_iiwa/joint_states"); //to know the robot status
//for(int i=0; i<7; i++ ) q_(i) = msg_r->position[i];        
//for(int i=0; i<7; i++ ) dq_(i) = msg_r->velocity[i];

//out = KUKA_INVKIN::forw_dy2(q_,dq_,v_); //predict dynamic for a inter-steps out = [q_ dq_]
//for (int j=0;j<7;j++) q_(j) = out(j); //unpack position
//for (int j=0;j<7;j++) dq_(j) = out(j+7); //unpack velocity
cout<<"q: "<<q_.transpose()<<endl;
cout<<"dq: "<<dq_.transpose()<<endl;

for(int i=0; i<7; i++ ) cmd_pub_q.data[i] = q_(i);
for(int i=0; i<7; i++ ) cmd_pub_v.data[i] = v_(i);
_my_pub_q.publish(cmd_pub_q);
_my_pub_v.publish(cmd_pub_v);

  */