#!groovy
//
// Jenkinsfile to control automated testing using https://jenkins.tjhei.info
//
// To run the tests locally, you can execute the following in a burnman working directory:
//
// docker run --rm -v "$(pwd):/src" tjhei/burnman:v4 /bin/bash -c "cd /src && ./test.sh"


pipeline
{
  agent none

  stages
  {

    stage("check")
    {
      agent
      {
        docker
        {
          image 'dealii/indent'
        }
      }

      post { cleanup { cleanWs() } }

      stages
      {
        stage("init")
        {
          steps
          {
            sh 'echo ${BRANCH_NAME} - PR #${CHANGE_ID}'
            githubNotify context: 'CI', description: 'need ready to test label and /rebuild',  status: 'PENDING'
          }
        }

        stage("permission")
        {
          // skip this if it is not a pull request:
          when { expression { env.CHANGE_ID != null } }
          steps
          {
            // For /rebuild to work you need to:
            // 1) select "issue comment" to be delivered in the github webhook setting
            // 2) install "GitHub PR Comment Build Plugin" on Jenkins
            // 3) in project settings select "add property" "Trigger build on pr comment" with
            //    the phrase ".*/rebuild.*" (without quotes)
            sh '''
               wget -q -O - https://api.github.com/repos/geodynamics/burnman/issues/${CHANGE_ID}/labels | grep 'ready to test' || \
               { echo "This commit will only be tested when it has the label 'ready to test'. Trigger a rebuild by adding a comment that contains '/rebuild'..."; exit 1; }
            '''
          }
        }

        stage("set status")
        {
          steps
          {
            githubNotify context: 'CI', description: 'running tests...',  status: 'PENDING'
          }
        }


        stage("indent")
        {
          steps
          {
            sh '''
               echo "indenting: TODO"
            '''
            sh 'git diff > changes.diff'
            archiveArtifacts artifacts: 'changes.diff', fingerprint: true
            sh '''
               git diff --exit-code || \
               { echo "Please check indentation, see artifacts in the top right corner!"; exit 1; }
            '''
            githubNotify context: 'indent', description: '',  status: 'SUCCESS'
          }

          post { failure {
          githubNotify context: 'indent', description: '',  status: 'FAILURE'
          }}

        }
      }
    }

    stage('Test')
    {
      agent
      {
        docker
        {
          image 'tjhei/burnman:v4'
        }
      }
      post { cleanup { cleanWs() } }


      stages
      {

        stage('python 2')
        {
          steps
          {
            timeout(time: 1, unit: 'HOURS')
            {
            sh '''#!/bin/bash
               PYTHON=python ./test.sh
            '''
            }
          }
        }

        stage('python 3')
        {
          steps
          {
            timeout(time: 1, unit: 'HOURS')
            {
            sh '''#!/bin/bash
               PYTHON=python3 ./test.sh
            '''
            }
          }
        }
      }
    }


    stage("finalize")
    {
      agent none

      steps
      {
        githubNotify context: 'CI', description: 'OK',  status: 'SUCCESS'
      }
    }
  }
}
