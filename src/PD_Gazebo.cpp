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
#include <Eigen/Dense>
#include <ros/package.h>
#include <tf/tf.h>
#include <tf_conversions/tf_eigen.h>

#include <kdl_parser/kdl_parser.hpp>
#include <kdl/chainfksolvervel_recursive.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolvervel_pinv.hpp>
#include <kdl/chainfksolverpos_recursive.hpp>
#include <kdl/chainiksolverpos_nr.hpp>
#include <kdl/chainjnttojacsolver.hpp>
#include <kdl/chainjnttojacdotsolver.hpp>
#include <kdl/chaindynparam.hpp>

#include "my_roscpp_library/my_roscpp_library.h"
#include "my_roscpp_function/my_roscpp_function.h"

#include <chrono>

using namespace Eigen;
using namespace std;
using namespace std::chrono;


class KUKA_INVDYN {
	public:
		KUKA_INVDYN();
		void run();
		bool init_robot_model();
		void joint_states_cb( sensor_msgs::JointState );
		void ctrl_loop();
                void ref_Callback(const std_msgs::Float64MultiArray::ConstPtr& msg);
 


	private:
		ros::NodeHandle _nh;
		KDL::Tree iiwa_tree;


		KDL::Chain _k_chain;
                VectorXd ref_q = VectorXd::Zero(7);


	
		ros::Subscriber _js_sub;
                ros::Subscriber _v_ERG;

		KDL::JntArray *_q_in;
		KDL::JntArray *_dq_in;
		bool _first_js;

		ros::Publisher _cmd_pub[7];

		KDL::ChainDynParam *_dyn_param;
};


bool KUKA_INVDYN::init_robot_model() {
	std::string robot_desc_string;
	_nh.param("robot_description", robot_desc_string, std::string());
	if (!kdl_parser::treeFromString(robot_desc_string, iiwa_tree)){
		ROS_ERROR("Failed to construct kdl tree");
		return false;
	}

	std::string base_link = "lbr_iiwa_link_0";
	std::string tip_link  = "lbr_iiwa_link_ee";
	if ( !iiwa_tree.getChain(base_link, tip_link, _k_chain) ) return false;


	_q_in = new KDL::JntArray( _k_chain.getNrOfJoints() );
	_dq_in = new KDL::JntArray( _k_chain.getNrOfJoints() );
	
	_dyn_param = new KDL::ChainDynParam(_k_chain,KDL::Vector(0,0,-9.81));


	return true;
}


KUKA_INVDYN::KUKA_INVDYN() {

	if (!init_robot_model()) exit(1); 
	ROS_INFO("Robot tree correctly loaded from parameter server!");

	cout << "Joints and segments: " << iiwa_tree.getNrOfJoints() << " - " << iiwa_tree.getNrOfSegments() << endl;
 
	_js_sub = _nh.subscribe("/lbr_iiwa/joint_states", 0, &KUKA_INVDYN::joint_states_cb, this);
	_v_ERG = _nh.subscribe("/ref", 0, &KUKA_INVDYN::ref_Callback, this);
	
	_cmd_pub[0] = _nh.advertise< std_msgs::Float64 > ("/lbr_iiwa/lbr_iiwa_joint_1_effort_controller/command", 0);
	_cmd_pub[1] = _nh.advertise< std_msgs::Float64 > ("/lbr_iiwa/lbr_iiwa_joint_2_effort_controller/command", 0);
	_cmd_pub[2] = _nh.advertise< std_msgs::Float64 > ("/lbr_iiwa/lbr_iiwa_joint_3_effort_controller/command", 0);
	_cmd_pub[3] = _nh.advertise< std_msgs::Float64 > ("/lbr_iiwa/lbr_iiwa_joint_4_effort_controller/command", 0);
	_cmd_pub[4] = _nh.advertise< std_msgs::Float64 > ("/lbr_iiwa/lbr_iiwa_joint_5_effort_controller/command", 0);
	_cmd_pub[5] = _nh.advertise< std_msgs::Float64 > ("/lbr_iiwa/lbr_iiwa_joint_6_effort_controller/command", 0);
	_cmd_pub[6] = _nh.advertise< std_msgs::Float64 > ("/lbr_iiwa/lbr_iiwa_joint_7_effort_controller/command", 0);

	_first_js = false;
	
}


void KUKA_INVDYN::ref_Callback(const std_msgs::Float64MultiArray::ConstPtr& msg) {

  for(int i = 0; i<7; i++) ref_q(i)= msg->data[i];
   cout<<"ref: "<<ref_q.transpose()<<endl;     
}


void KUKA_INVDYN::joint_states_cb( sensor_msgs::JointState js ) {

	for(int i=0; i<7; i++ ) { 
		_q_in->data[i] = js.position[i];
		_dq_in->data[i] = js.velocity[i];
		
	}
	_first_js = true;

}





void KUKA_INVDYN::ctrl_loop() {


ros::Rate r(1000);
KDL::JntArray grav_(7);

VectorXd q_des = VectorXd::Zero(7);


VectorXd _q = VectorXd::Zero(7);
MatrixXd Kp = MatrixXd::Zero(7,7);
MatrixXd Kd = MatrixXd::Zero(7,7);
VectorXd tao = VectorXd::Zero(7);
VectorXd q_tilda = VectorXd::Zero(7);
VectorXd _dq = VectorXd::Zero(7);


std_msgs::Float64 cmd[7];


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


while(!_first_js) usleep(0.1);
//for(int i = 0; i<7; i++) q_des(i) = _q_in->data[i];
for(int i = 0; i<7; i++) ref_q(i) = _q_in->data[i];
cout<<"first joint received: "<<q_des.transpose()<<endl;


while( ros::ok() ) { 



for(int i = 0; i<7; i++) _q(i) = _q_in->data[i] ;
for(int i = 0; i<7; i++)  _dq(i) = _dq_in->data[i];

q_tilda = ref_q - _q;
_dyn_param->JntToGravity(*_q_in, grav_);


tao = Kp*q_tilda - Kd*_dq + grav_.data;
for(int i=0; i<7; i++ )cmd[i].data = tao(i);
for(int i=0; i<7; i++ ) _cmd_pub[i].publish( cmd[i] );


r.sleep();
}

}


void KUKA_INVDYN::run() {


	boost::thread ctrl_loop_t ( &KUKA_INVDYN::ctrl_loop, this);
	ros::spin();	

}




int main(int argc, char** argv) {
        
        cout<<"StartComp"<<endl;
	ros::init(argc, argv, "kdl_node");
	KUKA_INVDYN ik;
	ik.run();

	return 0;
}
