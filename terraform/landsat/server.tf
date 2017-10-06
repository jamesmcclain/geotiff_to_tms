variable "region" {
  type        = "string"
  description = "AWS Region"
  default     = "us-west-2"
}

variable "ami" {
  type        = "string"
  description = "AMI (e.g. ami-e689729e)"
  default     = "ami-e689729e"
}

variable "key_name" {
  type        = "string"
  description = "The name of the EC2 secret key (primarily for SSH access)"
}

provider "aws" {
  region = "${var.region}"
}

resource "aws_instance" "landsat" {
  ami             = "${var.ami}"
  instance_type   = "t2.micro"
  key_name        = "${var.key_name}"
  security_groups = ["${aws_security_group.landsat.name}"]

  tags {
    Name = "LandSat"
  }
}

resource "aws_security_group" "landsat" {

  ingress {
    from_port   = "22"
    to_port     = "22"
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }

  ingress {
    from_port   = "8001"
    to_port     = "8001"
    protocol    = "tcp"
    cidr_blocks = ["0.0.0.0/0"]
  }

  egress {
    from_port   = 0
    to_port     = 0
    protocol    = "-1"
    cidr_blocks = ["0.0.0.0/0"]
  }

  lifecycle {
    create_before_destroy = true
  }
}

output "emr-master" {
  value = "${aws_instance.landsat.public_dns}"
}
