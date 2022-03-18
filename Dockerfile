ARG ROS_DISTRO=melodic

FROM ros:$ROS_DISTRO
LABEL description="geographiclib_conversions"
SHELL ["/bin/bash", "-c"]
WORKDIR /catkin_ws/src/geographiclib_conversions


# 1. Install basic requirements
RUN apt-get update                                                              &&  \
    apt-get install -y  ros-$ROS_DISTRO-catkin                                      \
                        python3-pip                                                 \
                        python3-catkin-tools
RUN if [[ "$ROS_DISTRO" = "melodic" ]] ; then apt-get install -y python-pip python-catkin-tools ; fi

# 2. Copy repository
COPY ./ ./

# 3. Setup packages
RUN ./scripts/install.sh

# 4. Build ROS
RUN source /opt/ros/$ROS_DISTRO/setup.bash && cd ../../ && catkin build

CMD echo "main process has been started"                                        &&  \
    source /opt/ros/$ROS_DISTRO/setup.bash                                      &&  \
    source /catkin_ws/devel/setup.bash                                          &&  \
    echo "container has been finished"
