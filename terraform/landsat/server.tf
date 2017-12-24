variable "region" {
  type        = "string"
  description = "AWS Region"
  default     = "us-west-2"
}

variable "ami" {
  type        = "string"
  description = "AMI (e.g. ami-bf4193c7)"
  default     = "ami-bf4193c7"
}

variable "key_name" {
  type        = "string"
  description = "The name of the EC2 secret key (primarily for SSH access)"
}

variable "instance_type" {
  type        = "string"
  description = "The type of instance (e.g. t2.micro)"
  default     = "t2.micro"
}

variable "tms_cidr" {
  type        = "string"
  description = "The CIDR block from which web connections are allowed"
  default     = "0.0.0.0/0"
}

provider "aws" {
  region = "${var.region}"
}

resource "aws_security_group" "tms" {

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
    cidr_blocks = ["${var.tms_cidr}"]
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

resource "aws_instance" "tms" {
  ami             = "${var.ami}"
  instance_type   = "${var.instance_type}"
  key_name        = "${var.key_name}"
  security_groups = ["${aws_security_group.tms.name}"]

  tags {
    Name = "TMS"
  }
}


output "emr-master" {
  value = "${aws_instance.tms.public_dns}"
}
